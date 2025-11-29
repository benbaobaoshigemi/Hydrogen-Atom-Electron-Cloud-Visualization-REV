/**
 * 手势控制模块 - 基于手部位置追踪的直接映射系统
 * VERSION 2.8 - 纯惯性控制系统
 * 
 * 交互逻辑（优化后的设计）：
 * ==========================================
 * 
 * ✊ 单手握拳拖拽 = 旋转视角
 *    - 握拳后，手部移动直接映射为相机旋转
 *    - 左右移动 = 水平旋转，上下移动 = 垂直旋转
 *    - 就像用手捏着物体旋转一样直观
 * 
 * 🖐️ 张开手掌 = 无动作
 *    - 张开手掌本身不触发任何动作
 *    - 旋转停止完全由惯性自然衰减决定
 * 
 * ✊/🖐️ 握拳/松手 = 状态切换
 *    - 握拳和松手动作本身不代表任何操作
 *    - 只是开始/结束追踪手部移动
 * 
 * ⏳ 惯性运行中 = 锁定操作
 *    - 当视角还在惯性旋转时，不识别任何新操作
 *    - 必须等待旋转完全停止后才能开始新的控制
 * 
 * 🤏🤏 双手捏合 = 缩放（需要两手分开）
 *    - 双手都捏合且相距足够远时，靠近/远离控制缩放
 *    - 惯性运行中时缩放也被锁定
 * 
 * ==========================================
 */

import {
    FilesetResolver,
    HandLandmarker
} from "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.0";

// ============ 版本信息 ============
console.log('%c[gesture.js] v3.0 - 纯惯性控制 + 颜色状态系统', 'color: lime; font-size: 14px; font-weight: bold;');
// ==================================

// ========================================
// 配置参数（从 gesture_config.json 加载）
// ========================================
let CONFIG = {
    // 捏合检测
    pinchStartThreshold: 0.055,
    pinchEndThreshold: 0.10,
    
    // 旋转控制
    rotationSensitivity: 2.5,
    deadzone: 0.004,
    maxDelta: 0.12,
    
    // 平滑滤波
    smoothingFactor: 0.35,
    
    // 复位保护
    resetProtection: true,
    directionChangeThreshold: 0.75,
    resetCooldownFrames: 8,
    velocityDecayOnReset: 0.1,
    
    // 释放缓冲
    releaseBufferFrames: 4,
    
    // 惯性系统
    friction: 0.97,
    minVelocity: 0.005,  // 大幅提高停止阈值，肉眼看上去停了就算停了
    inertiaBoost: 1.5,
    
    // 双手缩放
    zoomSensitivity: 2.5,
    minHandsSeparation: 0.25,
    handsOverlapThreshold: 0.15
};

// 加载配置文件
async function loadConfig() {
    try {
        const response = await fetch('gesture_config.json?t=' + Date.now());
        if (response.ok) {
            const json = await response.json();
            // 解析配置（支持带 description 的格式）
            for (const key in json) {
                if (key.startsWith('_') || key.startsWith('=')) continue;
                const val = json[key];
                if (val && typeof val === 'object' && 'value' in val) {
                    CONFIG[key] = val.value;
                } else if (typeof val !== 'object') {
                    CONFIG[key] = val;
                }
            }
            console.log('%c[gesture.js] 配置文件加载成功', 'color: cyan;', CONFIG);
        }
    } catch (e) {
        console.warn('[gesture.js] 配置文件加载失败，使用默认值:', e.message);
    }
}

// 立即加载配置
loadConfig();

// ========================================
// 模块状态
// ========================================
let handLandmarker = undefined;
let runningMode = "VIDEO";
let webcamRunning = false;
const video = document.createElement("video");
video.autoplay = true;
video.playsInline = true;
video.style.display = "none";
document.body.appendChild(video);

let lastVideoTime = -1;
let results = undefined;
let isModelLoading = false;
let modelLoadError = null;

// ========================================
// 交互状态
// ========================================

// 单手旋转状态
let isDragging = false;
let lastHandPosition = null;  // { x, y } 上一帧手部位置（原始值，不平滑）
let dragHandIndex = -1;       // 正在拖拽的手的索引

// 双手缩放状态
let isPinchZooming = false;
let lastPinchDistance = null;

// 惯性系统
let rotationVelocity = { x: 0, y: 0 };

// ========================================
// 速度平滑系统 - 对速度而非位置进行平滑
// ========================================
let smoothedVelocity = { x: 0, y: 0 };  // 平滑后的速度
let lastRawPosition = null;              // 上一帧原始位置（用于计算速度）

// ========================================
// 释放缓冲系统 - 防止松开时误操作
// ========================================
let releaseBufferFrames = 0;     // 释放缓冲帧计数

// ========================================
// 握拳状态滞后 - 避免在边界抖动
// ========================================
let wasFist = false;             // 上一帧是否握拳

// ========================================
// 复位保护系统 - 防止手复位时拽回视角
// ========================================
let lastMoveDirection = null;    // 上一帧移动方向 { x, y }
let resetCooldown = 0;           // 复位冷却帧计数

// ========================================
// Canvas 用于绘制手部骨架和状态
// ========================================
const canvasElement = document.createElement("canvas");
canvasElement.id = "gesture-canvas";
canvasElement.style.position = "absolute";
canvasElement.style.top = "0";
canvasElement.style.left = "0";
canvasElement.style.width = "100%";
canvasElement.style.height = "100%";
canvasElement.style.zIndex = "9998";
canvasElement.style.pointerEvents = "none";
document.body.appendChild(canvasElement);
const canvasCtx = canvasElement.getContext("2d");

// ========================================
// 状态提示函数
// 颜色系统：绿色=可操控，橙色=惯性滑动中，红色=拖拽中
// ========================================
function updateStatus(message, state = 'waiting') {
    const popup = document.getElementById('gesture-status-popup');
    const text = document.getElementById('gesture-status-text');
    const indicator = document.querySelector('.status-indicator');
    
    if (popup && text) {
        popup.style.display = 'flex';
        text.innerText = message;
        if (indicator) {
            switch(state) {
                case 'ready':     indicator.style.backgroundColor = '#00ff00'; break;  // 绿色 - 可操控
                case 'inertia':   indicator.style.backgroundColor = '#ff8800'; break;  // 橙色 - 惯性滑动中
                case 'dragging':  indicator.style.backgroundColor = '#ff0000'; break;  // 红色 - 拖拽中
                case 'waiting':   indicator.style.backgroundColor = 'yellow'; break;   // 黄色 - 等待
                case 'error':     indicator.style.backgroundColor = 'red'; break;
                default:          indicator.style.backgroundColor = 'yellow';
            }
        }
    }
    console.log(`[Gesture] ${message}`);
}

// ========================================
// 使用 MediaPipe HandLandmarker（纯手部追踪，更稳定）
// ========================================
const createHandLandmarker = async () => {
    if (isModelLoading || handLandmarker) return;
    isModelLoading = true;
    modelLoadError = null;
    
    try {
        console.log("正在加载手部追踪模型 (HandLandmarker)...");
        const btn = document.getElementById('gesture-control-btn');
        if(btn) btn.title = "正在加载模型...";

        const vision = await FilesetResolver.forVisionTasks(
            "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.0/wasm"
        );
        
        handLandmarker = await HandLandmarker.createFromOptions(vision, {
            baseOptions: {
                modelAssetPath: "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task",
                delegate: "GPU"
            },
            runningMode: runningMode,
            numHands: 2
        });
        
        console.log("HandLandmarker 模型加载完成");
        console.log("交互说明: 捏合拖拽旋转 | 张开手掌释放 | 双手捏合缩放");
        isModelLoading = false;
        
        if(btn) btn.title = "手势控制";
        
    } catch (error) {
        console.error("模型加载失败:", error);
        modelLoadError = error;
        isModelLoading = false;
        updateStatus("模型加载失败，请检查网络", 'error');
    }
};

createHandLandmarker();

// ========================================
// 公开 API
// ========================================
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Gesture = {};

window.ElectronCloud.Gesture.start = async function() {
    if (window.location.protocol === 'file:') {
        alert("错误：手势识别无法在本地文件模式(file://)下运行。\n\n请运行文件夹中的 'start_server.py' 脚本，或使用 VS Code Live Server 插件。");
        return;
    }

    if (modelLoadError) {
        alert("模型加载失败，请检查网络连接后刷新页面重试。\n" + modelLoadError.message);
        createHandLandmarker();
        return;
    }

    if (!handLandmarker) {
        if (isModelLoading) {
            updateStatus("正在初始化模型，请稍候...", 'waiting');
            const checkLoad = setInterval(() => {
                if (handLandmarker) {
                    clearInterval(checkLoad);
                    window.ElectronCloud.Gesture.start();
                } else if (modelLoadError) {
                    clearInterval(checkLoad);
                    updateStatus("模型加载失败", 'error');
                }
            }, 500);
            return;
        } else {
            await createHandLandmarker();
            if (!handLandmarker) return;
        }
    }
    
    if (webcamRunning === true) {
        return;
    }

    const constraints = {
        video: {
            width: 640,
            height: 480
        }
    };

    try {
        const stream = await navigator.mediaDevices.getUserMedia(constraints);
        video.srcObject = stream;
        video.onloadedmetadata = () => {
            video.play();
        };
        video.addEventListener("loadeddata", predictWebcam);
        webcamRunning = true;
        
        const popup = document.getElementById('gesture-status-popup');
        if(popup) {
            popup.style.display = 'flex';
            popup.style.zIndex = '9999';
        }
        updateStatus("🖐️ 张开 = 待机\n🤏 捏合拖拽 = 旋转", 'ready');
        
        requestAnimationFrame(physicsLoop);
        
    } catch (err) {
        console.error("Error accessing webcam:", err);
        let msg = "无法访问摄像头。";
        if (err.name === 'NotAllowedError') {
            msg += "请允许浏览器访问摄像头权限。";
        } else if (err.name === 'NotFoundError') {
            msg += "未检测到摄像头设备。";
        } else if (window.location.protocol !== 'https:' && window.location.hostname !== 'localhost' && window.location.hostname !== '127.0.0.1') {
            msg += "浏览器限制：摄像头只能在 HTTPS 或 localhost 下使用。";
        }
        alert(msg);
        updateStatus("摄像头启动失败", 'error');
    }
};

window.ElectronCloud.Gesture.stop = function() {
    webcamRunning = false;
    rotationVelocity = { x: 0, y: 0 };
    isDragging = false;
    isPinchZooming = false;
    lastHandPosition = null;
    lastPinchDistance = null;
    wasFist = false;
    lastRawPosition = null;
    smoothedVelocity = { x: 0, y: 0 };
    
    if(video.srcObject) {
        const tracks = video.srcObject.getTracks();
        tracks.forEach(track => track.stop());
        video.srcObject = null;
    }
    
    const popup = document.getElementById('gesture-status-popup');
    if(popup) popup.style.display = 'none';
    
    // 清除 canvas
    if (canvasElement && canvasCtx) {
        canvasCtx.clearRect(0, 0, canvasElement.width, canvasElement.height);
    }
    
    // 更新按钮状态
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.classList.remove('gesture-active');
        btn.title = '手势控制';
    }
    
    console.log('[Gesture] 手势控制已停止');
};

// 检查手势识别是否正在运行
window.ElectronCloud.Gesture.isRunning = function() {
    return webcamRunning;
};

// ========================================
// 检查惯性是否正在运行
// ========================================
function isInertiaRunning() {
    return Math.abs(rotationVelocity.x) > CONFIG.minVelocity || Math.abs(rotationVelocity.y) > CONFIG.minVelocity;
}

// ========================================
// 物理循环：应用惯性旋转
// ========================================
let inertiaLogCounter = 0;
function physicsLoop() {
    if (!webcamRunning) return;

    // 只有在非拖拽状态下才应用惯性
    if (!isDragging && !isPinchZooming) {
        if (isInertiaRunning()) {
            // 每20帧打印一次
            if (inertiaLogCounter++ % 20 === 0) {
                console.log('%c[INERTIA]', 'color: gold;', 
                    'velocity:', rotationVelocity.x.toFixed(5), rotationVelocity.y.toFixed(5),
                    'friction:', CONFIG.friction);
            }
            applyQuaternionRotation(rotationVelocity.x, rotationVelocity.y);
            
            // 应用阻尼
            rotationVelocity.x *= CONFIG.friction;
            rotationVelocity.y *= CONFIG.friction;
            
            if (Math.abs(rotationVelocity.x) < CONFIG.minVelocity) rotationVelocity.x = 0;
            if (Math.abs(rotationVelocity.y) < CONFIG.minVelocity) rotationVelocity.y = 0;
        }
    }

    requestAnimationFrame(physicsLoop);
}

// ========================================
// 四元数旋转实现 - 使用屏幕/窗口参考系
// deltaX: 屏幕上的左右旋转 (绕相机的上方向向量)
// deltaY: 屏幕上的上下旋转 (绕相机的右方向向量)
// ========================================
function applyQuaternionRotation(deltaX, deltaY) {
    const state = window.ElectronCloud.state;
    if (!state || !state.camera || !state.controls) return;

    const camera = state.camera;
    const controls = state.controls;
    const target = controls.target;

    // 获取相机到目标的偏移向量
    const offset = new THREE.Vector3().subVectors(camera.position, target);
    
    // 使用屏幕/窗口参考系：
    // 左右旋转：绕相机的 "up" 向量 (屏幕的上方向)
    // 上下旋转：绕相机的 "right" 向量 (屏幕的右方向)
    
    // 获取相机的本地坐标轴
    const cameraUp = camera.up.clone().normalize();
    const cameraRight = new THREE.Vector3();
    cameraRight.crossVectors(camera.getWorldDirection(new THREE.Vector3()), cameraUp).normalize();
    
    // 创建旋转四元数
    // 绕相机 up 轴旋转 (屏幕上的左右)
    const qHorizontal = new THREE.Quaternion().setFromAxisAngle(cameraUp, -deltaX);
    
    // 绕相机 right 轴旋转 (屏幕上的上下)
    const qVertical = new THREE.Quaternion().setFromAxisAngle(cameraRight, -deltaY);
    
    // 组合旋转
    const quaternion = new THREE.Quaternion();
    quaternion.multiplyQuaternions(qVertical, qHorizontal);
    
    // 应用旋转到偏移向量
    offset.applyQuaternion(quaternion);
    
    // 更新相机位置
    camera.position.copy(target).add(offset);
    
    // 更新相机的上向量以保持正确的朝向
    camera.up.applyQuaternion(quaternion);
    camera.up.normalize();
    
    // 让相机看向目标
    camera.lookAt(target);
    
    // 同步 OrbitControls 状态
    controls.update();
}

// ========================================
// 预测循环
// ========================================
async function predictWebcam() {
    if (!webcamRunning) return;
    
    try {
        if (video.paused || video.ended) {
            await video.play();
        }

        let startTimeMs = performance.now();
        if (lastVideoTime !== video.currentTime) {
            lastVideoTime = video.currentTime;
            results = handLandmarker.detectForVideo(video, startTimeMs);
        }
        
        processHandTracking(results);
    } catch (error) {
        console.error("手势识别循环错误:", error);
    }
    
    if (webcamRunning) {
        window.requestAnimationFrame(predictWebcam);
    }
}

// ========================================
// 判断是否捏合（拇指和食指靠近）- 带滞后
// ========================================
// ========================================
// 速度计算与平滑 - 直接使用原始位置，对速度进行平滑
// ========================================
function calculateSmoothedVelocity(currentPos) {
    if (!lastRawPosition) {
        lastRawPosition = { x: currentPos.x, y: currentPos.y };
        smoothedVelocity = { x: 0, y: 0 };
        console.log('%c[VEL INIT]', 'color: magenta;', 
            'pos:', currentPos.x.toFixed(4), currentPos.y.toFixed(4));
        return { x: 0, y: 0 };
    }
    
    // 计算原始速度（位置差）
    const rawVelX = currentPos.x - lastRawPosition.x;
    const rawVelY = currentPos.y - lastRawPosition.y;
    
    // 对速度进行 EMA 平滑
    const alpha = CONFIG.smoothingFactor;
    smoothedVelocity.x = alpha * rawVelX + (1 - alpha) * smoothedVelocity.x;
    smoothedVelocity.y = alpha * rawVelY + (1 - alpha) * smoothedVelocity.y;
    
    // 更新上一帧位置
    lastRawPosition.x = currentPos.x;
    lastRawPosition.y = currentPos.y;
    
    // 每10帧打印一次
    if (Math.random() < 0.1) {
        console.log('%c[VEL]', 'color: magenta;', 
            'raw:', rawVelX.toFixed(5), rawVelY.toFixed(5),
            '→ smooth:', smoothedVelocity.x.toFixed(5), smoothedVelocity.y.toFixed(5));
    }
    
    return { x: smoothedVelocity.x, y: smoothedVelocity.y };
}

// ========================================
// 握拳检测（带滞后）- 检测所有手指是否弯曲
// ========================================
let fistLogCounter = 0;
function isFistWithHysteresis(landmarks, wasFistBefore) {
    // 手指关键点索引:
    // 拇指: 1-4 (CMC, MCP, IP, TIP)
    // 食指: 5-8 (MCP, PIP, DIP, TIP)
    // 中指: 9-12
    // 无名指: 13-16
    // 小指: 17-20
    // 手腕: 0
    
    const wrist = landmarks[0];
    const palmBase = landmarks[9];  // 中指MCP作为手掌基准
    
    // 计算每个手指尖到手腕的距离，与手指根部到手腕的距离比较
    // 如果指尖比指根更近手腕，说明手指弯曲了
    
    // 食指
    const indexTip = landmarks[8];
    const indexMCP = landmarks[5];
    const indexBent = distance2D(indexTip, wrist) < distance2D(indexMCP, wrist) * 1.1;
    
    // 中指
    const middleTip = landmarks[12];
    const middleMCP = landmarks[9];
    const middleBent = distance2D(middleTip, wrist) < distance2D(middleMCP, wrist) * 1.1;
    
    // 无名指
    const ringTip = landmarks[16];
    const ringMCP = landmarks[13];
    const ringBent = distance2D(ringTip, wrist) < distance2D(ringMCP, wrist) * 1.1;
    
    // 小指
    const pinkyTip = landmarks[20];
    const pinkyMCP = landmarks[17];
    const pinkyBent = distance2D(pinkyTip, wrist) < distance2D(pinkyMCP, wrist) * 1.1;
    
    // 拇指 - 检查是否收拢（靠近手掌）
    const thumbTip = landmarks[4];
    const thumbBent = distance2D(thumbTip, palmBase) < 0.15;
    
    // 计算弯曲的手指数量
    const bentCount = [indexBent, middleBent, ringBent, pinkyBent].filter(b => b).length;
    
    // 滞后判定
    let result;
    if (wasFistBefore) {
        // 如果之前是握拳，需要至少2个手指伸直才算松开
        result = bentCount >= 2;
    } else {
        // 如果之前不是握拳，需要至少3个手指弯曲才算握拳
        result = bentCount >= 3;
    }
    
    // 每15帧打印一次，或者状态变化时
    const shouldLog = (fistLogCounter++ % 15 === 0) || (result !== wasFistBefore);
    if (shouldLog) {
        console.log('%c[FIST]', result ? 'color: green; font-weight: bold;' : 'color: gray;',
            'bent:', bentCount, '/4',
            'thumb:', thumbBent,
            'result:', result,
            'wasF:', wasFistBefore);
    }
    
    return result;
}

// 辅助函数：计算2D距离
function distance2D(p1, p2) {
    return Math.sqrt(Math.pow(p1.x - p2.x, 2) + Math.pow(p1.y - p2.y, 2));
}

// ========================================
// 简单捏合检测（用于双手缩放）
// ========================================
function isPinching(landmarks) {
    const thumbTip = landmarks[4];
    const indexTip = landmarks[8];
    
    const distance = Math.sqrt(
        Math.pow(thumbTip.x - indexTip.x, 2) + 
        Math.pow(thumbTip.y - indexTip.y, 2)
    );
    
    return distance < 0.07;
}

// ========================================
// 获取手部边界框（用于检测重叠）
// ========================================
function getHandBoundingBox(landmarks) {
    let minX = 1, maxX = 0, minY = 1, maxY = 0;
    for (const lm of landmarks) {
        minX = Math.min(minX, lm.x);
        maxX = Math.max(maxX, lm.x);
        minY = Math.min(minY, lm.y);
        maxY = Math.max(maxY, lm.y);
    }
    return { minX, maxX, minY, maxY };
}

// ========================================
// 检查两只手是否分离（不重叠）
// ========================================
function areHandsSeparated(hands) {
    if (hands.length < 2) return false;
    
    const center1 = getHandCenter(hands[0]);
    const center2 = getHandCenter(hands[1]);
    
    // 计算两手中心距离
    const centerDistance = Math.sqrt(
        Math.pow(center1.x - center2.x, 2) + 
        Math.pow(center1.y - center2.y, 2)
    );
    
    // 两手中心必须距离足够远
    if (centerDistance < CONFIG.minHandsSeparation) {
        return false;
    }
    
    // 检查边界框是否重叠
    const box1 = getHandBoundingBox(hands[0]);
    const box2 = getHandBoundingBox(hands[1]);
    
    // 计算重叠区域
    const overlapX = Math.max(0, Math.min(box1.maxX, box2.maxX) - Math.max(box1.minX, box2.minX));
    const overlapY = Math.max(0, Math.min(box1.maxY, box2.maxY) - Math.max(box1.minY, box2.minY));
    const overlapArea = overlapX * overlapY;
    
    // 计算较小手的面积
    const area1 = (box1.maxX - box1.minX) * (box1.maxY - box1.minY);
    const area2 = (box2.maxX - box2.minX) * (box2.maxY - box2.minY);
    const smallerArea = Math.min(area1, area2);
    
    // 如果重叠面积超过较小手面积的阈值，认为重叠
    if (smallerArea > 0 && overlapArea / smallerArea > CONFIG.handsOverlapThreshold) {
        return false;
    }
    
    return true;
}

// ========================================
// 获取手部中心位置（使用手腕和中指MCP的中点）
// ========================================
function getHandCenter(landmarks) {
    const wrist = landmarks[0];
    const middleMCP = landmarks[9];
    
    return {
        x: (wrist.x + middleMCP.x) / 2,
        y: (wrist.y + middleMCP.y) / 2
    };
}

// ========================================
// 核心处理逻辑
// ========================================
function processHandTracking(results) {
    const popupText = document.getElementById('gesture-status-text');
    
    // 检查惯性是否正在运行
    const inertiaActive = !isDragging && !isPinchZooming && isInertiaRunning();
    
    // 无手检测
    if (!results || !results.landmarks || results.landmarks.length === 0) {
        // 根据状态显示不同颜色
        if (isDragging) {
            updateStatus("✊ 拖拽中...", 'dragging');
        } else if (isInertiaRunning()) {
            updateStatus("⏳ 惯性旋转中...", 'inertia');
        } else {
            updateStatus("等待检测手部...", 'waiting');
        }
        
        if (canvasElement && canvasCtx) {
            canvasCtx.clearRect(0, 0, canvasElement.width, canvasElement.height);
        }
        
        // 停止所有交互，但保留惯性（不清零 rotationVelocity）
        if (isDragging) {
            isDragging = false;
            lastRawPosition = null;
            smoothedVelocity = { x: 0, y: 0 };
            wasFist = false;
        }
        isPinchZooming = false;
        lastPinchDistance = null;
        return;
    }

    const hands = results.landmarks;
    const handedness = results.handedness;
    
    // 绘制手部骨架（始终绘制）
    drawHandSkeleton(hands, handedness);
    
    // ========================================
    // 惯性运行中时，如果握拳就直接开始控制（清除惯性）
    // 只有张开手掌时才继续惯性滑行
    // ========================================
    if (inertiaActive) {
        // 识别当前手势状态
        const landmarks = hands[0];
        const isFist = isFistWithHysteresis(landmarks, wasFist);
        wasFist = isFist;
        
        if (isFist) {
            // 握拳了！直接清除惯性并开始控制
            rotationVelocity = { x: 0, y: 0 };  // 清除惯性
            // 不 return，继续执行下面的握拳拖拽逻辑
        } else {
            // 张开手掌，继续惯性滑行
            updateStatus("⏳ 惯性旋转中...\n🖐️ 张开手掌", 'inertia');
            return;
        }
    }
    
    // ========================================
    // 双手捏合缩放（需要两只手分开且都捏合）
    // ========================================
    if (hands.length >= 2 && areHandsSeparated(hands)) {
        const hand1Pinching = isPinching(hands[0]);
        const hand2Pinching = isPinching(hands[1]);
        
        if (hand1Pinching && hand2Pinching) {
            const center1 = getHandCenter(hands[0]);
            const center2 = getHandCenter(hands[1]);
            
            const currentDistance = Math.sqrt(
                Math.pow(center1.x - center2.x, 2) + 
                Math.pow(center1.y - center2.y, 2)
            );
            
            if (lastPinchDistance !== null) {
                const delta = currentDistance - lastPinchDistance;
                if (Math.abs(delta) > 0.005) {
                    handleZoom(delta * CONFIG.zoomSensitivity);
                    updateStatus(delta > 0 ? "🤏🤏 放大中..." : "🤏🤏 缩小中...", 'active');
                }
            } else {
                updateStatus("🤏🤏 双手缩放模式", 'active');
            }
            
            lastPinchDistance = currentDistance;
            isPinchZooming = true;
            isDragging = false;
            lastRawPosition = null;
            smoothedVelocity = { x: 0, y: 0 };
            return;
        }
    }
    
    // 重置双手缩放状态
    if (isPinchZooming) {
        isPinchZooming = false;
        lastPinchDistance = null;
    }
    
    // ========================================
    // 单手握拳控制旋转
    // ========================================
    const landmarks = hands[0];
    const rawHandCenter = getHandCenter(landmarks);
    const isFist = isFistWithHysteresis(landmarks, wasFist);
    
    // 释放缓冲：刚松开时忽略几帧
    if (releaseBufferFrames > 0) {
        releaseBufferFrames--;
        updateStatus("🖐️ 释放中...", 'ready');
        wasFist = false;
        lastMoveDirection = null;
        return;
    }
    
    // 复位冷却中
    if (resetCooldown > 0) {
        resetCooldown--;
        // 冷却期间仍然更新位置，但不应用旋转
        if (isFist) {
            // 更新原始位置追踪，但不计算速度
            lastRawPosition = { x: rawHandCenter.x, y: rawHandCenter.y };
            smoothedVelocity = { x: 0, y: 0 };
            lastMoveDirection = null;  // 重置方向
        }
        updateStatus("✊ 复位中...", 'active');
        return;
    }
    
    if (isFist) {
        // 握拳状态 = 开始或继续拖拽
        wasFist = true;
        
        if (!isDragging) {
            // 开始拖拽 - 初始化速度系统
            isDragging = true;
            lastRawPosition = { x: rawHandCenter.x, y: rawHandCenter.y };
            smoothedVelocity = { x: 0, y: 0 };
            lastMoveDirection = null;
            rotationVelocity = { x: 0, y: 0 }; // 清除惯性
            console.log('%c[DRAG START]', 'color: yellow; font-weight: bold;', 
                'pos:', rawHandCenter.x.toFixed(4), rawHandCenter.y.toFixed(4));
            // 红色 - 拖拽中
            updateStatus("✊ 握拳拖拽中...\n移动手部旋转视角", 'dragging');
        } else {
            // 继续拖拽 - 使用速度平滑系统
            // 计算平滑后的速度
            const velocity = calculateSmoothedVelocity(rawHandCenter);
            const velMag = Math.sqrt(velocity.x * velocity.x + velocity.y * velocity.y);
            
            // ========================================
            // 复位保护：检测方向突变
            // ========================================
            if (CONFIG.resetProtection && lastMoveDirection && velMag > 0.0005) {
                const lastMag = Math.sqrt(lastMoveDirection.x * lastMoveDirection.x + lastMoveDirection.y * lastMoveDirection.y);
                if (lastMag > 0.0005) {
                    const dot = (velocity.x * lastMoveDirection.x + velocity.y * lastMoveDirection.y) / (velMag * lastMag);
                    
                    console.log('%c[DIR CHECK]', 'color: orange;', 
                        'dot:', dot.toFixed(3), 
                        'threshold:', CONFIG.directionChangeThreshold,
                        'velMag:', velMag.toFixed(5));
                    
                    // 如果方向变化超过阈值（dot < threshold），认为是复位动作
                    if (dot < CONFIG.directionChangeThreshold) {
                        // 触发复位保护
                        console.log('%c[RESET TRIGGERED!]', 'color: red; font-weight: bold; font-size: 14px;',
                            'cooldown:', CONFIG.resetCooldownFrames);
                        resetCooldown = CONFIG.resetCooldownFrames;
                        rotationVelocity.x *= CONFIG.velocityDecayOnReset;
                        rotationVelocity.y *= CONFIG.velocityDecayOnReset;
                        // 重置速度系统
                        smoothedVelocity = { x: 0, y: 0 };
                        lastMoveDirection = null;
                        updateStatus("✊ 检测到复位...", 'dragging');
                        return;
                    }
                }
            }
            
            // 更新移动方向（只有速度足够大才更新）
            if (velMag > 0.0005) {
                lastMoveDirection = { x: velocity.x, y: velocity.y };
            }
            
            // 应用灵敏度（镜像处理：摄像头左右是反的）
            const deltaX = -velocity.x * CONFIG.rotationSensitivity;
            const deltaY = velocity.y * CONFIG.rotationSensitivity;
            
            // 死区过滤
            if (Math.abs(deltaX) > CONFIG.deadzone || Math.abs(deltaY) > CONFIG.deadzone) {
                // 限幅
                const clampedDeltaX = Math.max(-CONFIG.maxDelta, Math.min(CONFIG.maxDelta, deltaX));
                const clampedDeltaY = Math.max(-CONFIG.maxDelta, Math.min(CONFIG.maxDelta, deltaY));
                
                console.log('%c[ROTATE]', 'color: lime;', 
                    'vel:', velocity.x.toFixed(5), velocity.y.toFixed(5),
                    'delta:', clampedDeltaX.toFixed(4), clampedDeltaY.toFixed(4));
                
                applyQuaternionRotation(clampedDeltaX, clampedDeltaY);
                
                // 记录速度用于惯性（应用放大系数）
                const boost = CONFIG.inertiaBoost || 1.0;
                rotationVelocity.x = clampedDeltaX * boost;
                rotationVelocity.y = clampedDeltaY * boost;
                
                // 红色 - 拖拽中
                updateStatus("🔴 拖拽中...\n移动手部旋转视角", 'dragging');
            } else {
                // 在死区内时，让速度平滑器更快衰减
                smoothedVelocity.x *= 0.5;
                smoothedVelocity.y *= 0.5;
                // 红色 - 拖拽中（静止）
                updateStatus("🔴 握拳中...\n移动手部旋转视角", 'dragging');
            }
        }
    } else {
        // 手掌张开 = 释放控制（但不意味着停止）
        if (isDragging) {
            console.log('%c[DRAG END]', 'color: cyan; font-weight: bold;',
                'velocity:', rotationVelocity.x.toFixed(4), rotationVelocity.y.toFixed(4));
            isDragging = false;
            lastRawPosition = null;  // 重置速度系统
            smoothedVelocity = { x: 0, y: 0 };
            lastMoveDirection = null;
            releaseBufferFrames = CONFIG.releaseBufferFrames;  // 启动释放缓冲
            // 保留惯性速度 - 旋转停止只由惯性自然衰减决定
        }
        
        wasFist = false;
        // 张开手掌不触发任何动作，惯性会继续运行
        // 颜色：橙色=惯性滑动中，绿色=可操控
        if (isInertiaRunning()) {
            updateStatus("⏳ 惯性旋转中...\n🖐️ 等待停止", 'inertia');
        } else {
            updateStatus("🟢 就绪\n握拳开始控制", 'ready');
        }
    }
}

// ========================================
// 绘制手部骨架
// ========================================
function drawHandSkeleton(landmarksArray, handednessArray) {
    if (!canvasElement || !canvasCtx) return;

    const displayWidth = canvasElement.clientWidth;
    const displayHeight = canvasElement.clientHeight;

    if (canvasElement.width !== displayWidth || canvasElement.height !== displayHeight) {
        canvasElement.width = displayWidth;
        canvasElement.height = displayHeight;
    }
    
    canvasCtx.clearRect(0, 0, canvasElement.width, canvasElement.height);
    
    if (!landmarksArray || landmarksArray.length === 0) return;

    const width = canvasElement.width;
    const height = canvasElement.height;

    for (let handIndex = 0; handIndex < landmarksArray.length; handIndex++) {
        const landmarks = landmarksArray[handIndex];
        const isPinch = isPinching(landmarks);
        
        // 根据手势状态选择颜色
        let strokeColor;
        let statusEmoji;
        if (isPinch) {
            strokeColor = "#00FF00";  // 绿色 - 捏合（控制中）
            statusEmoji = "🤏";
        } else {
            strokeColor = "#FFFFFF";  // 白色 - 张开（待机）
            statusEmoji = "🖐️";
        }
        
        canvasCtx.lineWidth = 3;
        canvasCtx.strokeStyle = strokeColor;

        // 手指连接
        const connections = [
            [0, 1], [1, 2], [2, 3], [3, 4],       // 拇指
            [0, 5], [5, 6], [6, 7], [7, 8],       // 食指
            [0, 9], [9, 10], [10, 11], [11, 12],  // 中指
            [0, 13], [13, 14], [14, 15], [15, 16], // 无名指
            [0, 17], [17, 18], [18, 19], [19, 20], // 小指
            [5, 9], [9, 13], [13, 17]             // 手掌
        ];

        for (const [start, end] of connections) {
            const p1 = landmarks[start];
            const p2 = landmarks[end];
            
            // 镜像显示
            const x1 = (1 - p1.x) * width;
            const y1 = p1.y * height;
            const x2 = (1 - p2.x) * width;
            const y2 = p2.y * height;
            
            canvasCtx.beginPath();
            canvasCtx.moveTo(x1, y1);
            canvasCtx.lineTo(x2, y2);
            canvasCtx.stroke();
        }

        // 绘制关键点
        canvasCtx.fillStyle = strokeColor;
        for (const landmark of landmarks) {
            const x = (1 - landmark.x) * width;
            const y = landmark.y * height;
            
            canvasCtx.beginPath();
            canvasCtx.arc(x, y, 4, 0, 2 * Math.PI);
            canvasCtx.fill();
        }
        
        // 显示手势状态
        const wrist = landmarks[0];
        const labelX = (1 - wrist.x) * width;
        const labelY = wrist.y * height - 20;
        
        canvasCtx.fillStyle = strokeColor;
        canvasCtx.font = "bold 24px Arial";
        canvasCtx.textAlign = "center";
        canvasCtx.fillText(statusEmoji, labelX, labelY);
    }
}

// ========================================
// 缩放处理
// ========================================
function handleZoom(delta) {
    const state = window.ElectronCloud.state;
    if (!state || !state.controls || !state.camera) {
        return;
    }
    
    const controls = state.controls;
    const camera = state.camera;

    const offset = new THREE.Vector3().copy(camera.position).sub(controls.target);
    const currentDist = offset.length();
    let newDist = currentDist;

    if (delta > 0) {
        newDist = currentDist / (1 + Math.abs(delta));
    } else {
        newDist = currentDist * (1 + Math.abs(delta));
    }

    newDist = Math.max(1, Math.min(newDist, 500));

    offset.setLength(newDist);
    camera.position.copy(controls.target).add(offset);
    
    controls.update();
}
