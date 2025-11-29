/**
 * gesture_v2.js - 重构版手势控制系统
 * 
 * 版本: 4.0 - 捏合旋转模式
 * 
 * 交互设计（参考 GitHub 上成功案例）:
 * - 🤏 单手捏合 + 移动 = 旋转视角（可无限次抬起放下实现360°+旋转）
 * - 🤏🤏 双手捏合 + 拉开/靠近 = 缩放
 * - 🖐️ 张开手掌 = 释放（惯性继续）
 * 
 * 核心改进:
 * 1. 捏合检测比握拳检测稳定得多（只需判断拇指-食指距离）
 * 2. 简化状态机：只有 idle, rotating, zooming 三种状态
 * 3. 移除复杂的方向突变检测、复位保护等补丁
 */

// ========================================
// 命名空间初始化
// ========================================
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Gesture = window.ElectronCloud.Gesture || {};

// ========================================
// 配置参数
// ========================================
const CONFIG = {
    // 捏合检测（相对距离 = 手指距离/手部大小）
    pinchStartThreshold: 0.35,    // 开始捏合（相对距离，约35%手掌宽度）
    pinchEndThreshold: 0.18,      // 结束捏合（缩小到1/3，约18%手掌宽度）
    pinchReleaseVelocity: 0.08,   // 手指分开速度阈值（快速分开=意图松开）
    
    // 旋转控制
    rotationSensitivity: 5.0,     // 旋转灵敏度
    deadzone: 0.001,              // 死区（更灵敏）
    maxDelta: 0.15,               // 单帧最大位移限幅
    
    // 平滑
    smoothingFactor: 0.35,        // EMA 平滑系数 (降低，更平滑)
    pinchDistanceSmoothing: 0.5,  // 捏合距离平滑系数
    
    // 惯性
    friction: 0.98,               // 阻尼系数
    minVelocity: 0.0001,          // 惯性停止阈值
    inertiaBoost: 1.2,            // 惯性初始速度放大
    
    // 缩放
    zoomSensitivity: 3.0,         // 缩放灵敏度
    minHandSeparation: 0.15,      // 双手最小间距（防止误检）
    
    // 释放缓冲
    releaseBufferFrames: 3,       // 松开后忽略的帧数
};

// ========================================
// 状态变量
// ========================================
let handLandmarker = null;
let webcamRunning = false;
let video = null;
let canvasElement = null;
let canvasCtx = null;
let lastVideoTime = -1;
let results = null;

// 简化的状态机
const STATE = {
    IDLE: 'idle',           // 空闲/待机
    ROTATING: 'rotating',   // 单手捏合旋转中
    ZOOMING: 'zooming'      // 双手捏合缩放中
};
let currentState = STATE.IDLE;

// 捏合状态（带滞后 + 意图识别）
let isPinching = false;          // 当前是否捏合
let lastPinchPosition = null;    // 上一帧捏合位置
let smoothedPosition = null;     // 平滑后的位置
let lastPinchDist = null;        // 上一帧捏合距离（用于意图识别）
let smoothedPinchDist = null;    // 平滑后的捏合距离

// 惯性
let rotationVelocity = { x: 0, y: 0 };

// 双手缩放
let lastPinchDistance = null;

// 释放缓冲
let releaseBuffer = 0;

// ========================================
// 动态创建 UI 元素
// ========================================
function createGestureElements() {
    // 创建视频元素（用于摄像头输入）
    if (!document.getElementById('gesture-video')) {
        video = document.createElement('video');
        video.id = 'gesture-video';
        video.style.display = 'none';
        video.setAttribute('playsinline', '');
        video.setAttribute('autoplay', '');
        document.body.appendChild(video);
    } else {
        video = document.getElementById('gesture-video');
    }
    
    // 创建 canvas 元素（用于绘制手部骨架）
    if (!document.getElementById('gesture-canvas')) {
        canvasElement = document.createElement('canvas');
        canvasElement.id = 'gesture-canvas';
        canvasElement.style.position = 'absolute';
        canvasElement.style.top = '0';
        canvasElement.style.left = '0';
        canvasElement.style.width = '100%';
        canvasElement.style.height = '100%';
        canvasElement.style.zIndex = '9998';
        canvasElement.style.pointerEvents = 'none';
        document.body.appendChild(canvasElement);
    } else {
        canvasElement = document.getElementById('gesture-canvas');
    }
    
    canvasCtx = canvasElement.getContext('2d');
    
    console.log('[Gesture v4] UI 元素已创建');
}

// ========================================
// 初始化
// ========================================
async function initializeHandLandmarker() {
    const { HandLandmarker, FilesetResolver } = await import(
        "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.0"
    );

    const vision = await FilesetResolver.forVisionTasks(
        "https://cdn.jsdelivr.net/npm/@mediapipe/tasks-vision@0.10.0/wasm"
    );

    handLandmarker = await HandLandmarker.createFromOptions(vision, {
        baseOptions: {
            modelAssetPath: "https://storage.googleapis.com/mediapipe-models/hand_landmarker/hand_landmarker/float16/1/hand_landmarker.task",
            delegate: "GPU"
        },
        runningMode: "VIDEO",
        numHands: 2
    });

    console.log('[Gesture v4] HandLandmarker 初始化完成');
}

// ========================================
// UI 状态更新
// ========================================
function updateStatus(text, state = 'ready') {
    const popup = document.getElementById('gesture-status-popup');
    const textEl = document.getElementById('gesture-status-text');
    
    if (textEl) {
        textEl.textContent = text;
    }
    
    if (popup) {
        popup.className = 'gesture-status-popup';
        popup.classList.add(`status-${state}`);
    }
}

// ========================================
// 计算手部大小（用于归一化距离）
// 使用手腕到中指MCP的距离作为参考
// ========================================
function getHandSize(landmarks) {
    const wrist = landmarks[0];
    const middleMCP = landmarks[9];
    
    return Math.hypot(
        middleMCP.x - wrist.x,
        middleMCP.y - wrist.y
    );
}

// ========================================
// 捏合检测（带滞后 + 意图识别 + 相对距离）
// ========================================
function checkPinch(landmarks, wasPinching) {
    const thumbTip = landmarks[4];
    const indexTip = landmarks[8];
    
    // 计算绝对距离
    const absoluteDistance = Math.hypot(
        thumbTip.x - indexTip.x,
        thumbTip.y - indexTip.y
    );
    
    // 计算手部大小并归一化距离
    const handSize = getHandSize(landmarks);
    const distance = handSize > 0.01 ? absoluteDistance / handSize : absoluteDistance;
    
    // 平滑捏合距离
    if (smoothedPinchDist === null) {
        smoothedPinchDist = distance;
    } else {
        smoothedPinchDist = CONFIG.pinchDistanceSmoothing * distance + 
                            (1 - CONFIG.pinchDistanceSmoothing) * smoothedPinchDist;
    }
    
    // 意图识别：检测手指分开速度
    let releaseIntent = false;
    if (lastPinchDist !== null && wasPinching) {
        const distVelocity = smoothedPinchDist - lastPinchDist;
        // 手指快速分开 = 意图松开
        if (distVelocity > CONFIG.pinchReleaseVelocity) {
            releaseIntent = true;
            console.log('[Gesture] 检测到松开意图, 速度:', distVelocity.toFixed(4));
        }
    }
    lastPinchDist = smoothedPinchDist;
    
    // 调试日志（偶尔输出）
    if (Math.random() < 0.05) {
        console.log('[Pinch] 相对距离:', distance.toFixed(3), 
                    '手大小:', handSize.toFixed(3),
                    '阈值:', CONFIG.pinchStartThreshold);
    }
    
    // 滞后判定 + 意图识别
    if (wasPinching) {
        // 之前在捏合：距离超过阈值 或 检测到松开意图 都算松开
        if (releaseIntent || distance > CONFIG.pinchEndThreshold) {
            return false;
        }
        return true;
    } else {
        // 之前未捏合，需要距离小于开始阈值才算捏合
        return distance < CONFIG.pinchStartThreshold;
    }
}

// ========================================
// 获取捏合点位置（拇指和食指的中点）
// ========================================
function getPinchPosition(landmarks) {
    const thumbTip = landmarks[4];
    const indexTip = landmarks[8];
    
    return {
        x: (thumbTip.x + indexTip.x) / 2,
        y: (thumbTip.y + indexTip.y) / 2
    };
}

// ========================================
// 获取手部中心位置
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
// 计算两手距离（用于缩放）
// ========================================
function getHandsDistance(hand1, hand2) {
    const center1 = getHandCenter(hand1);
    const center2 = getHandCenter(hand2);
    
    return Math.hypot(
        center1.x - center2.x,
        center1.y - center2.y
    );
}

// ========================================
// 检查两只手是否足够分离
// ========================================
function areHandsSeparated(hand1, hand2) {
    return getHandsDistance(hand1, hand2) > CONFIG.minHandSeparation;
}

// ========================================
// 四元数旋转实现
// ========================================
function applyRotation(deltaX, deltaY) {
    const state = window.ElectronCloud.state;
    if (!state || !state.camera || !state.controls) return;

    const camera = state.camera;
    const controls = state.controls;
    const target = controls.target;

    const offset = new THREE.Vector3().subVectors(camera.position, target);
    
    // 获取相机坐标系
    const cameraUp = camera.up.clone().normalize();
    const cameraRight = new THREE.Vector3();
    cameraRight.crossVectors(camera.getWorldDirection(new THREE.Vector3()), cameraUp).normalize();
    
    // 创建旋转四元数
    const qHorizontal = new THREE.Quaternion().setFromAxisAngle(cameraUp, -deltaX);
    const qVertical = new THREE.Quaternion().setFromAxisAngle(cameraRight, -deltaY);
    
    const quaternion = new THREE.Quaternion();
    quaternion.multiplyQuaternions(qVertical, qHorizontal);
    
    // 应用旋转
    offset.applyQuaternion(quaternion);
    camera.position.copy(target).add(offset);
    camera.up.applyQuaternion(quaternion).normalize();
    camera.lookAt(target);
    controls.update();
}

// ========================================
// 缩放处理
// ========================================
function applyZoom(delta) {
    const state = window.ElectronCloud.state;
    if (!state || !state.controls || !state.camera) return;
    
    const controls = state.controls;
    const camera = state.camera;

    const offset = new THREE.Vector3().copy(camera.position).sub(controls.target);
    const currentDist = offset.length();
    
    let newDist;
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

// ========================================
// 惯性是否运行中
// ========================================
function isInertiaActive() {
    return Math.abs(rotationVelocity.x) > CONFIG.minVelocity || 
           Math.abs(rotationVelocity.y) > CONFIG.minVelocity;
}

// ========================================
// 物理循环（惯性）
// ========================================
function physicsLoop() {
    if (!webcamRunning) return;

    // 只有在空闲状态才应用惯性
    if (currentState === STATE.IDLE && isInertiaActive()) {
        applyRotation(rotationVelocity.x, rotationVelocity.y);
        
        // 阻尼衰减
        rotationVelocity.x *= CONFIG.friction;
        rotationVelocity.y *= CONFIG.friction;
        
        if (Math.abs(rotationVelocity.x) < CONFIG.minVelocity) rotationVelocity.x = 0;
        if (Math.abs(rotationVelocity.y) < CONFIG.minVelocity) rotationVelocity.y = 0;
    }

    requestAnimationFrame(physicsLoop);
}

// ========================================
// 核心处理逻辑
// ========================================
function processHands(results) {
    // 无手检测
    if (!results || !results.landmarks || results.landmarks.length === 0) {
        // 如果正在操作，停止并保留惯性
        if (currentState === STATE.ROTATING) {
            currentState = STATE.IDLE;
            lastPinchPosition = null;
            smoothedPosition = null;
            releaseBuffer = CONFIG.releaseBufferFrames;
        }
        if (currentState === STATE.ZOOMING) {
            currentState = STATE.IDLE;
            lastPinchDistance = null;
        }
        
        if (canvasElement && canvasCtx) {
            canvasCtx.clearRect(0, 0, canvasElement.width, canvasElement.height);
        }
        
        if (isInertiaActive()) {
            updateStatus("⏳ 惯性旋转中...", 'inertia');
        } else {
            updateStatus("等待检测手部...", 'waiting');
        }
        return;
    }

    const hands = results.landmarks;
    
    // 绘制手部
    drawHands(hands, results.handedness);
    
    // 释放缓冲期
    if (releaseBuffer > 0) {
        releaseBuffer--;
        return;
    }
    
    // ========================================
    // 优先处理双手缩放
    // ========================================
    if (hands.length >= 2) {
        const hand1 = hands[0];
        const hand2 = hands[1];
        
        const hand1Pinching = checkPinch(hand1, currentState === STATE.ZOOMING);
        const hand2Pinching = checkPinch(hand2, currentState === STATE.ZOOMING);
        
        if (hand1Pinching && hand2Pinching && areHandsSeparated(hand1, hand2)) {
            const currentDistance = getHandsDistance(hand1, hand2);
            
            if (currentState === STATE.ZOOMING && lastPinchDistance !== null) {
                const delta = currentDistance - lastPinchDistance;
                if (Math.abs(delta) > 0.005) {
                    applyZoom(delta * CONFIG.zoomSensitivity);
                    updateStatus(delta > 0 ? "🤏🤏 放大中..." : "🤏🤏 缩小中...", 'active');
                }
            } else {
                updateStatus("🤏🤏 双手缩放模式", 'active');
            }
            
            lastPinchDistance = currentDistance;
            currentState = STATE.ZOOMING;
            
            // 重置单手状态
            isPinching = false;
            lastPinchPosition = null;
            smoothedPosition = null;
            rotationVelocity = { x: 0, y: 0 };
            return;
        }
    }
    
    // 退出缩放状态
    if (currentState === STATE.ZOOMING) {
        currentState = STATE.IDLE;
        lastPinchDistance = null;
    }
    
    // ========================================
    // 单手捏合旋转
    // ========================================
    const hand = hands[0];
    const nowPinching = checkPinch(hand, isPinching);
    
    if (nowPinching) {
        const pinchPos = getPinchPosition(hand);
        
        if (!isPinching) {
            // 刚开始捏合 - 初始化
            isPinching = true;
            currentState = STATE.ROTATING;
            lastPinchPosition = pinchPos;
            smoothedPosition = pinchPos;
            rotationVelocity = { x: 0, y: 0 };  // 清除惯性
            
            console.log('[Gesture] 开始捏合旋转');
            updateStatus("🤏 捏合旋转中...\n移动手部旋转视角", 'active');
        } else {
            // 继续捏合 - 计算移动
            
            // EMA 平滑当前位置
            const smoothX = CONFIG.smoothingFactor * pinchPos.x + 
                           (1 - CONFIG.smoothingFactor) * smoothedPosition.x;
            const smoothY = CONFIG.smoothingFactor * pinchPos.y + 
                           (1 - CONFIG.smoothingFactor) * smoothedPosition.y;
            
            // 计算位移（相对于上一帧平滑位置）
            const deltaX = smoothX - smoothedPosition.x;
            const deltaY = smoothY - smoothedPosition.y;
            
            // 更新平滑位置
            smoothedPosition.x = smoothX;
            smoothedPosition.y = smoothY;
            
            // 计算位移大小
            const deltaMag = Math.sqrt(deltaX * deltaX + deltaY * deltaY);
            
            // 调试日志（每10帧输出一次）
            if (Math.random() < 0.1) {
                console.log('[Gesture] delta:', deltaX.toFixed(5), deltaY.toFixed(5), 
                            'mag:', deltaMag.toFixed(5), 'deadzone:', CONFIG.deadzone);
            }
            
            // 死区过滤
            if (deltaMag > CONFIG.deadzone) {
                // 限幅
                const clampedX = Math.max(-CONFIG.maxDelta, Math.min(CONFIG.maxDelta, deltaX));
                const clampedY = Math.max(-CONFIG.maxDelta, Math.min(CONFIG.maxDelta, deltaY));
                
                // 应用旋转（镜像：摄像头是镜像的）
                const rotX = -clampedX * CONFIG.rotationSensitivity;
                const rotY = clampedY * CONFIG.rotationSensitivity;
                
                console.log('[Gesture] 旋转:', rotX.toFixed(4), rotY.toFixed(4));
                
                applyRotation(rotX, rotY);
                
                // 记录速度用于惯性
                rotationVelocity.x = rotX * CONFIG.inertiaBoost;
                rotationVelocity.y = rotY * CONFIG.inertiaBoost;
                
                updateStatus("🤏 旋转中...", 'active');
            } else {
                updateStatus("🤏 捏合中（静止）", 'active');
            }
        }
    } else {
        // 松开捏合
        if (isPinching) {
            console.log('[Gesture] 松开捏合，惯性:', 
                rotationVelocity.x.toFixed(4), rotationVelocity.y.toFixed(4));
            
            isPinching = false;
            currentState = STATE.IDLE;
            lastPinchPosition = null;
            smoothedPosition = null;
            lastPinchDist = null;        // 重置意图识别
            smoothedPinchDist = null;
            releaseBuffer = CONFIG.releaseBufferFrames;
            // 保留惯性速度
        }
        
        if (isInertiaActive()) {
            updateStatus("⏳ 惯性旋转中...\n🖐️ 张开等待", 'inertia');
        } else {
            updateStatus("🖐️ 就绪\n捏合拇指食指开始", 'ready');
        }
    }
}

// ========================================
// 绘制手部骨架
// ========================================
function drawHands(landmarksArray, handednessArray) {
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
        const pinching = checkPinch(landmarks, isPinching);
        
        // 颜色：捏合时绿色，张开时白色
        const strokeColor = pinching ? "#00FF00" : "#FFFFFF";
        const emoji = pinching ? "🤏" : "🖐️";
        
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
        
        // 高亮拇指和食指（用于捏合）
        const thumbTip = landmarks[4];
        const indexTip = landmarks[8];
        
        canvasCtx.fillStyle = "#FFD700";  // 金色
        canvasCtx.beginPath();
        canvasCtx.arc((1 - thumbTip.x) * width, thumbTip.y * height, 6, 0, 2 * Math.PI);
        canvasCtx.fill();
        canvasCtx.beginPath();
        canvasCtx.arc((1 - indexTip.x) * width, indexTip.y * height, 6, 0, 2 * Math.PI);
        canvasCtx.fill();
        
        // 捏合时画连线
        if (pinching) {
            canvasCtx.strokeStyle = "#FFD700";
            canvasCtx.lineWidth = 2;
            canvasCtx.beginPath();
            canvasCtx.moveTo((1 - thumbTip.x) * width, thumbTip.y * height);
            canvasCtx.lineTo((1 - indexTip.x) * width, indexTip.y * height);
            canvasCtx.stroke();
        }
        
        // 状态 emoji
        const wrist = landmarks[0];
        const labelX = (1 - wrist.x) * width;
        const labelY = wrist.y * height - 20;
        
        canvasCtx.fillStyle = strokeColor;
        canvasCtx.font = "bold 24px Arial";
        canvasCtx.textAlign = "center";
        canvasCtx.fillText(emoji, labelX, labelY);
    }
}

// ========================================
// 视频帧处理循环
// ========================================
async function predictWebcam() {
    if (!webcamRunning) return;
    
    try {
        if (video.paused || video.ended) {
            await video.play();
        }

        const startTimeMs = performance.now();
        if (lastVideoTime !== video.currentTime) {
            lastVideoTime = video.currentTime;
            results = handLandmarker.detectForVideo(video, startTimeMs);
        }
        
        processHands(results);
    } catch (error) {
        console.error("手势识别错误:", error);
    }
    
    if (webcamRunning) {
        requestAnimationFrame(predictWebcam);
    }
}

// ========================================
// 公开 API
// ========================================

/**
 * 启动手势控制
 */
window.ElectronCloud.Gesture.start = async function() {
    if (webcamRunning) {
        console.log('[Gesture] 已在运行中');
        return;
    }

    console.log('[Gesture v4] 启动手势控制（捏合旋转模式）');
    
    // 动态创建 UI 元素
    createGestureElements();
    
    if (!video || !canvasElement) {
        console.error('[Gesture] UI 元素创建失败');
        return;
    }
    
    // 初始化 HandLandmarker
    if (!handLandmarker) {
        updateStatus("加载手势模型...", 'waiting');
        await initializeHandLandmarker();
    }
    
    // 启动摄像头
    try {
        updateStatus("启动摄像头...", 'waiting');
        
        const stream = await navigator.mediaDevices.getUserMedia({
            video: {
                width: { ideal: 640 },
                height: { ideal: 480 },
                facingMode: "user"
            }
        });
        
        video.srcObject = stream;
        await video.play();
        
        webcamRunning = true;
        
        // 重置状态
        currentState = STATE.IDLE;
        isPinching = false;
        lastPinchPosition = null;
        smoothedPosition = null;
        rotationVelocity = { x: 0, y: 0 };
        lastPinchDistance = null;
        releaseBuffer = 0;
        
        // 启动循环
        requestAnimationFrame(predictWebcam);
        requestAnimationFrame(physicsLoop);
        
        // 更新 UI
        const btn = document.getElementById('gesture-control-btn');
        if (btn) {
            btn.classList.add('gesture-active');
            btn.title = '手势控制（运行中）';
        }
        
        const popup = document.getElementById('gesture-status-popup');
        if (popup) {
            popup.style.display = 'flex';
            popup.style.zIndex = '9999';
        }
        
        updateStatus("🖐️ 就绪\n捏合拇指食指开始", 'ready');
        
        console.log('[Gesture v4] 手势控制已启动');
        
    } catch (err) {
        console.error("摄像头启动失败:", err);
        let msg = "无法访问摄像头。";
        if (err.name === 'NotAllowedError') {
            msg += "请允许浏览器访问摄像头权限。";
        } else if (err.name === 'NotFoundError') {
            msg += "未检测到摄像头设备。";
        }
        alert(msg);
        updateStatus("摄像头启动失败", 'error');
    }
};

/**
 * 停止手势控制
 */
window.ElectronCloud.Gesture.stop = function() {
    webcamRunning = false;
    
    // 重置所有状态
    currentState = STATE.IDLE;
    isPinching = false;
    lastPinchPosition = null;
    smoothedPosition = null;
    rotationVelocity = { x: 0, y: 0 };
    lastPinchDistance = null;
    
    if (video && video.srcObject) {
        const tracks = video.srcObject.getTracks();
        tracks.forEach(track => track.stop());
        video.srcObject = null;
    }
    
    const popup = document.getElementById('gesture-status-popup');
    if (popup) popup.style.display = 'none';
    
    if (canvasElement && canvasCtx) {
        canvasCtx.clearRect(0, 0, canvasElement.width, canvasElement.height);
    }
    
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.classList.remove('gesture-active');
        btn.title = '手势控制';
    }
    
    console.log('[Gesture v4] 手势控制已停止');
};

/**
 * 检查是否正在运行
 */
window.ElectronCloud.Gesture.isRunning = function() {
    return webcamRunning;
};

console.log('[Gesture v4] 模块加载完成 - 捏合旋转模式');
