// Three.js 场景管理模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Scene = {};

// 初始化场景
window.ElectronCloud.Scene.init = function() {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    
    // 场景
    state.scene = new THREE.Scene();

    // 相机
    state.camera = new THREE.PerspectiveCamera(75, window.innerWidth / window.innerHeight, 0.1, 1000);
    state.camera.position.set(0, 0, 40);
    state.camera.lookAt(0, 0, 0);
    // 相机默认只渲染layer 0（点云）和 layer 1（坐标轴），但bloom只影响layer 0
    state.camera.layers.enable(0);
    state.camera.layers.enable(1);

    // 渲染器
    const container = document.getElementById('container');
    // 开启 alpha: true 以支持透明背景截图
    // 开启 preserveDrawingBuffer: true 以支持 toDataURL 导出图片
    state.renderer = new THREE.WebGLRenderer({ 
        antialias: true, 
        alpha: true,
        preserveDrawingBuffer: true 
    });
    
    // 【全屏画布方案】画布像素尺寸始终使用屏幕全屏尺寸
    // CSS 显示尺寸为窗口大小，这样导出时直接捕捉全屏画布
    const pixelRatio = Math.max(2, Math.min(window.devicePixelRatio, 3));
    state.fullscreenPixelRatio = pixelRatio; // 保存以便后续使用
    
    // 画布像素尺寸 = 屏幕全屏尺寸 * 像素比
    const canvasWidth = Math.floor(window.screen.width * pixelRatio);
    const canvasHeight = Math.floor(window.screen.height * pixelRatio);
    
    // 设置渲染器：第一个参数是CSS尺寸，setPixelRatio之后会乘以像素比
    // 为了得到全屏像素尺寸，我们先设置size为屏幕尺寸，再设置pixelRatio
    state.renderer.setSize(window.screen.width, window.screen.height);
    state.renderer.setPixelRatio(pixelRatio);
    
    // 然后通过CSS将显示尺寸限制为窗口大小
    state.renderer.domElement.style.width = '100%';
    state.renderer.domElement.style.height = '100%';
    
    console.log('渲染器像素比:', pixelRatio);
    console.log('画布像素尺寸 (全屏):', canvasWidth, 'x', canvasHeight);
    console.log('CSS显示尺寸 (窗口):', window.innerWidth, 'x', window.innerHeight);
    container.appendChild(state.renderer.domElement);

    // 初始化后期处理（Bloom辉光效果）
    window.ElectronCloud.Scene.initBloom();

    // 控制器
    state.controls = new THREE.OrbitControls(state.camera, state.renderer.domElement);
    state.controls.enableDamping = true;
    state.controls.dampingFactor = 0.05;
    // 禁用 OrbitControls 的内置旋转，使用自定义轨迹球旋转
    state.controls.enableRotate = false;
    
    // 初始化自定义轨迹球鼠标控制（完全自由旋转，无极点限制）
    window.ElectronCloud.Scene.initTrackballMouseControl();
    
    // 自动旋转状态
    state.autoRotate = {
        enabled: false,
        axis: new THREE.Vector3(0, 1, 0),
        speed: 1, // 默认速度为最大值的10%
        axisHelper: null
    };

    // 自定义坐标轴（默认不显示）
    state.customAxes = window.ElectronCloud.Scene.createCustomAxes(constants.AXES_BASE_SIZE);
    state.customAxes.visible = false;
    // 坐标轴放在layer 1，不受bloom影响
    state.customAxes.traverse(function(obj) {
        obj.layers.set(1);
    });
    state.scene.add(state.customAxes);
    
    // 权重模式状态（独立于闪烁模式）
    state.weightMode = false;
    
    // 闪烁模式状态
    state.heartbeat = {
        enabled: false,
        phase: 0,
        baseStrength: 0,
        maxBrightness: 200,
        mode: 'starry', // 'starry'(默认), 'breath', 'wave'(波浪)
        frequency: 50,     // 0-100，默认50%（中间值）
        starryPhases: null, // 星空模式的随机相位数组
        weightMode: false   // 兼容旧属性
    };

    // 监听窗口大小变化
    window.addEventListener('resize', window.ElectronCloud.Scene.onWindowResize, false);
    
    console.log('Three.js 场景初始化完成');
};

// 创建自定义坐标轴
window.ElectronCloud.Scene.createCustomAxes = function(size) {
    const axes = new THREE.Group();
    
    // 使用 Line2 实现粗线条并启用抗锯齿效果
    // 但如果 Line2 不可用，回退到 LineBasicMaterial
    const useThickLines = typeof THREE.LineMaterial !== 'undefined' && 
                          typeof THREE.LineGeometry !== 'undefined' &&
                          typeof THREE.Line2 !== 'undefined';
    
    // 基本白色材质（用于粗线或回退）
    const lineMaterial = new THREE.LineBasicMaterial({ 
        color: 0xffffff,
        linewidth: 2  // 注意：WebGL 中大多数浏览器忽略此属性
    });

    // 为每个轴创建独立的组，以便单独控制可见性
    const xAxisGroup = new THREE.Group();
    xAxisGroup.name = 'axis-x';
    const yAxisGroup = new THREE.Group();
    yAxisGroup.name = 'axis-y';
    const zAxisGroup = new THREE.Group();
    zAxisGroup.name = 'axis-z';

    // 轴线
    const pointsX = [new THREE.Vector3(-size, 0, 0), new THREE.Vector3(size, 0, 0)];
    const pointsY = [new THREE.Vector3(0, -size, 0), new THREE.Vector3(0, size, 0)];
    const pointsZ = [new THREE.Vector3(0, 0, -size), new THREE.Vector3(0, 0, size)];
    
    xAxisGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsX), lineMaterial));
    yAxisGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsY), lineMaterial));
    zAxisGroup.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsZ), lineMaterial));

    // 箭头 - 使用更多细分来减少锯齿
    const arrowGeom = new THREE.ConeGeometry(0.5, 2, 16); // 增加细分数从8到16
    const arrowMat = new THREE.MeshBasicMaterial({ color: 0xffffff });

    const arrowX = new THREE.Mesh(arrowGeom, arrowMat);
    arrowX.position.set(size, 0, 0);
    arrowX.rotation.z = -Math.PI / 2;
    xAxisGroup.add(arrowX);

    const arrowY = new THREE.Mesh(arrowGeom, arrowMat);
    arrowY.position.set(0, size, 0);
    yAxisGroup.add(arrowY);

    const arrowZ = new THREE.Mesh(arrowGeom, arrowMat);
    arrowZ.position.set(0, 0, size);
    arrowZ.rotation.x = Math.PI / 2;
    zAxisGroup.add(arrowZ);
    
    // 将各轴组添加到主组
    axes.add(xAxisGroup);
    axes.add(yAxisGroup);
    axes.add(zAxisGroup);
    
    // 标签
    const loader = new THREE.FontLoader();
    loader.load(
        'https://cdn.jsdelivr.net/npm/three@0.128.0/examples/fonts/helvetiker_regular.typeface.json',
        function (font) {
            const textMat = new THREE.MeshBasicMaterial({ color: 0xffffff });
            
            const textShapesX = font.generateShapes('X', 1.5);
            const textGeomX = new THREE.ShapeGeometry(textShapesX);
            const textMeshX = new THREE.Mesh(textGeomX, textMat);
            textMeshX.position.set(size + 2, 0, 0);
            xAxisGroup.add(textMeshX);

            const textShapesY = font.generateShapes('Y', 1.5);
            const textGeomY = new THREE.ShapeGeometry(textShapesY);
            const textMeshY = new THREE.Mesh(textGeomY, textMat);
            textMeshY.position.set(0, size + 2, 0);
            yAxisGroup.add(textMeshY);

            const textShapesZ = font.generateShapes('Z', 1.5);
            const textGeomZ = new THREE.ShapeGeometry(textShapesZ);
            const textMeshZ = new THREE.Mesh(textGeomZ, textMat);
            textMeshZ.position.set(0, 0, size + 2);
            zAxisGroup.add(textMeshZ);
        },
        undefined, // onProgress
        function (error) {
            console.warn('坐标轴字体加载失败，标签将不显示:', error);
        }
    );

    return axes;
};

// 设置某个轴的可见性
window.ElectronCloud.Scene.setAxisVisibility = function(axis, visible) {
    const state = window.ElectronCloud.state;
    if (!state.customAxes) return;
    
    const axisGroup = state.customAxes.getObjectByName('axis-' + axis);
    if (axisGroup) {
        axisGroup.visible = visible;
    }
};

// 重置所有轴的可见性
window.ElectronCloud.Scene.resetAxesVisibility = function() {
    const state = window.ElectronCloud.state;
    if (!state.customAxes) return;
    
    ['x', 'y', 'z'].forEach(axis => {
        const axisGroup = state.customAxes.getObjectByName('axis-' + axis);
        if (axisGroup) {
            axisGroup.visible = true;
        }
    });
};

// 处理窗口大小变化
window.ElectronCloud.Scene.onWindowResize = function() {
    const state = window.ElectronCloud.state;
    
    // 【全屏画布方案】只更新相机宽高比，不改变画布像素尺寸
    // 画布始终保持屏幕全屏尺寸，CSS自动缩放显示
    state.camera.aspect = window.innerWidth / window.innerHeight;
    state.camera.updateProjectionMatrix();
    
    // 不需要调用 renderer.setSize，因为 CSS style 已经设置为 100%
    // 画布像素尺寸保持不变（全屏尺寸）
    
    console.log('窗口调整，显示尺寸:', window.innerWidth, 'x', window.innerHeight,
                '画布像素尺寸保持:', state.renderer.domElement.width, 'x', state.renderer.domElement.height);
};

// 初始化Bloom辉光效果
window.ElectronCloud.Scene.initBloom = function() {
    const state = window.ElectronCloud.state;
    
    // 初始状态
    state.bloomEnabled = false;
    state.bloomStrength = 0;
    
    try {
        // 检查是否加载了后期处理库
        if (typeof THREE.EffectComposer === 'undefined' || 
            typeof THREE.RenderPass === 'undefined' ||
            typeof THREE.UnrealBloomPass === 'undefined') {
            console.warn('Bloom后期处理库未完全加载，辉光效果不可用');
            return;
        }
        
        // 【全屏画布方案】使用画布的实际像素尺寸（已经是全屏尺寸）
        const canvas = state.renderer.domElement;
        const width = canvas.width;
        const height = canvas.height;
        
        console.log('Bloom初始化使用全屏画布尺寸:', width, 'x', height);
        
        // 创建支持 alpha 通道的 RenderTarget（这是支持透明背景导出的关键）
        const renderTarget = new THREE.WebGLRenderTarget(width, height, {
            minFilter: THREE.LinearFilter,
            magFilter: THREE.LinearFilter,
            format: THREE.RGBAFormat,  // 使用 RGBA 格式以支持 alpha 通道
            stencilBuffer: false
        });
        
        // 创建EffectComposer，使用支持 alpha 的 renderTarget
        state.composer = new THREE.EffectComposer(state.renderer, renderTarget);
        
        // 渲染通道
        const renderPass = new THREE.RenderPass(state.scene, state.camera);
        renderPass.clearAlpha = 0;  // 清除时使用透明背景
        state.composer.addPass(renderPass);
        
        // 【关键】Bloom通道也使用实际像素分辨率
        state.bloomPass = new THREE.UnrealBloomPass(
            new THREE.Vector2(width, height),
            0,      // 初始强度 strength
            0.2,    // 半径 radius - 缩小以防止边缘被裁剪
            0       // 阈值 threshold - 设为0使所有颜色都能发光
        );
        state.composer.addPass(state.bloomPass);
        
        console.log('Bloom辉光效果初始化完成，分辨率:', width, 'x', height);
    } catch (err) {
        console.warn('Bloom辉光效果初始化失败:', err);
        state.composer = null;
        state.bloomPass = null;
    }
};

// 动画循环
window.ElectronCloud.Scene.animate = function() {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    
    state.animationFrameId = requestAnimationFrame(window.ElectronCloud.Scene.animate);

    // 如果正在绘制，则更新点
    if (state.isDrawing) {
        if (state.pointCount < state.MAX_POINTS) {
            // 优先使用 Web Worker 并行采样（如果可用）
            if (window.ElectronCloud.Sampling.isUsingWorkers && window.ElectronCloud.Sampling.isUsingWorkers()) {
                window.ElectronCloud.Sampling.submitSamplingTask();
            } else {
                // 降级到主线程采样
                window.ElectronCloud.Sampling.updatePoints();
            }
        }
        
        // 图表刷新：每秒刷新一次（1000ms），跳过第一秒
        // 【修复】将图表刷新移出 pointCount < MAX_POINTS 条件，确保采样全程都能刷新
        const now = performance.now();
        if (!state.lastChartUpdateTime) state.lastChartUpdateTime = 0;
        
        // 检查是否已过1秒延迟期（只跳过第一次刷新）
        const isDelayPassed = state.samplingStartTime && (now - state.samplingStartTime > 1000);
        
        // 【修复】放宽条件：只要有采样数据或点数增加，就尝试刷新图表
        // 这确保即使在高点数模式下，图表也能持续每秒刷新一次直到采样完成
        const hasData = state.radialSamples.length > 0 || state.pointCount > 0;
        
        if (isDelayPassed && hasData && (now - state.lastChartUpdateTime) >= 1000) {
            state.lastChartUpdateTime = now;
            
            // 只有当有径向数据时才更新图表数据
            if (state.radialSamples.length > 0) {
                window.ElectronCloud.Orbital.updateBackgroundChartData();
            }
            
            // 如果数据面板可见，实时刷新图表
            const dataPanel = document.getElementById('data-panel');
            if (dataPanel && !dataPanel.classList.contains('collapsed')) {
                window.ElectronCloud.Orbital.drawProbabilityChart(false);
            }
        }
        
        // 检查采样完成
        if (state.pointCount >= state.MAX_POINTS && !state.samplingCompleted) {
            window.ElectronCloud.Scene.onSamplingCompleted();
        }
    }

    // 确保坐标系状态和大小正确：实时更新坐标系
    if (state.customAxes && state.axesScaleFactor > 0) {
        // 使用公共函数计算缩放
        const result = window.ElectronCloud.UI.calculateAxesScale(state.axesScaleFactor);
        state.customAxes.visible = result.visible;
        
        if (result.visible) {
            // 只有当缩放比例发生变化时才更新（避免不必要的更新）
            const currentScale = state.customAxes.scale.x;
            if (Math.abs(currentScale - result.scale) > 0.001) {
                state.customAxes.scale.set(result.scale, result.scale, result.scale);
            }
        }
    } else if (state.customAxes) {
        // 比例系数为0时隐藏坐标系
        state.customAxes.visible = false;
    }
    
    // 自动旋转逻辑
    // 【关键修复】只有在点云存在时才执行自动旋转，避免采样前坐标轴独自旋转导致偏移
    if (state.autoRotate && state.autoRotate.enabled && state.autoRotate.speed > 0 && state.points) {
        // 在当前原子坐标系下旋转（局部坐标系）
        const axis = state.autoRotate.axis.clone().normalize();
        // 速度系数：最大速度10对应之前1.5的效果
        const angle = state.autoRotate.speed * 0.0015;
        
        // 绕指定轴旋转场景中的点云（使用局部坐标系）
        state.points.rotateOnAxis(axis, angle);
        
        // 角向叠加层也一起旋转
        if (state.angularOverlay) {
            state.angularOverlay.rotateOnAxis(axis, angle);
        }
        // 坐标轴也跟着旋转
        if (state.customAxes) {
            state.customAxes.rotateOnAxis(axis, angle);
        }
        
        // 累计旋转角度（用于录制功能）
        if (state.autoRotate.totalAngle === undefined) {
            state.autoRotate.totalAngle = 0;
        }
        state.autoRotate.totalAngle += angle;
        
        // 如果正在录制，显示进度
        if (state.isRecordingRotation) {
            const progress = (state.autoRotate.totalAngle / (Math.PI * 2) * 100).toFixed(1);
            // 已旋转一周（2π），自动停止录制
            if (state.autoRotate.totalAngle >= Math.PI * 2) {
                console.log('已旋转一周 (360°)，自动停止录制');
                state.isRecordingRotation = false; // 立即标记为停止，避免重复调用
                if (window.ElectronCloud.UI.stopRotationRecording) {
                    window.ElectronCloud.UI.stopRotationRecording();
                }
            }
        }
    }
    
    // 闪烁模式
    if (state.heartbeat && state.heartbeat.enabled && state.bloomPass) {
        const frequencyScale = 0.3 + (state.heartbeat.frequency / 100) * 1.4; // 0.3-1.7
        
        if (state.heartbeat.mode === 'starry') {
            // 星空模式：每个点独立随机闪烁，像星海一样
            state.heartbeat.phase += 0.06 * frequencyScale; // 原来0.02，现在翻3倍
            
            // 基础bloom强度 - 使用心跳模式的值的较高倍数，降低到原来的80%
            const maxStrength = (state.heartbeat.maxBrightness || 200) / 100 * 1.5 * 2.5 * 0.8;
            state.bloomPass.strength = maxStrength * 0.6; // 基础bloom
            
            // 更新每个点的亮度（通过临时修改颜色实现）
            if (state.points && state.points.geometry) {
                const colors = state.points.geometry.attributes.color;
                if (colors) {
                    const pointCount = state.pointCount || 0;
                    const time = state.heartbeat.phase;
                    
                    // 初始化随机相位数组（每个点一个独立的相位）
                    const maxPoints = state.MAX_POINTS || 100000;
                    if (!state.heartbeat.starryPhases || state.heartbeat.starryPhases.length < maxPoints) {
                        state.heartbeat.starryPhases = new Float32Array(maxPoints);
                        state.heartbeat.starryFrequencies = new Float32Array(maxPoints);
                        // 为每个点生成随机相位和频率
                        for (let i = 0; i < maxPoints; i++) {
                            state.heartbeat.starryPhases[i] = Math.random() * Math.PI * 2;
                            state.heartbeat.starryFrequencies[i] = 0.5 + Math.random() * 1.5;
                        }
                    }
                    
                    // 确保有基础颜色数组用于保存原始颜色
                    // 注意：新点的颜色在采样时已经同步到 baseColors 中
                    if (!state.baseColors || state.baseColors.length < maxPoints * 3) {
                        state.baseColors = new Float32Array(maxPoints * 3);
                        state.baseColorsCount = 0;
                    }
                    
                    // 如果星空模式刚刚启动，需要从当前颜色初始化已存在点的baseColors
                    // 但如果颜色已经被修改过，这里会读到错误的值
                    // 所以只在 baseColorsCount 为 0 时初始化
                    if (state.baseColorsCount === 0 && pointCount > 0) {
                        const colorArray = colors.array;
                        for (let i = 0; i < pointCount * 3; i++) {
                            state.baseColors[i] = colorArray[i];
                        }
                        state.baseColorsCount = pointCount;
                    }
                    
                    // 亮度范围：暗时允许接近透明（5%），亮时为心跳模式值的2.0倍（80%的2.5）
                    const minBrightness = 0.05;
                    const maxBrightnessMultiplier = 2.0;
                    
                    // 【权重模式】确保有概率密度权重数组
                    const positions = state.points.geometry.attributes.position;
                    if (state.heartbeat.weightMode && positions) {
                        if (!state.densityWeights || state.densityWeights.length < maxPoints) {
                            state.densityWeights = new Float32Array(maxPoints);
                            for (let i = 0; i < pointCount; i++) {
                                const i3 = i * 3;
                                const x = positions.array[i3];
                                const y = positions.array[i3 + 1];
                                const z = positions.array[i3 + 2];
                                const r = Math.sqrt(x * x + y * y + z * z);
                                state.densityWeights[i] = Math.exp(-r / 10) * 0.7 + 0.3;
                            }
                        }
                    }
                    
                    // 更新每个点的颜色亮度
                    const colorArray = colors.array;
                    
                    // 只处理已经有baseColors的点
                    const processCount = Math.min(pointCount, state.baseColorsCount);
                    
                    for (let i = 0; i < processCount; i++) {
                        const phase = state.heartbeat.starryPhases[i];
                        const freq = state.heartbeat.starryFrequencies[i];
                        
                        // 每个点独立的闪烁波形 (0-1)
                        const individualWave = 0.5 + 0.5 * Math.sin(time * freq + phase);
                        // 添加一些随机闪烁（偶尔变亮）
                        const sparkle = Math.random() < 0.002 ? 1.3 : 1.0;
                        
                        // 【权重模式】密集区域峰值亮度更高
                        let brightnessMultiplier = maxBrightnessMultiplier;
                        if (state.heartbeat.weightMode && state.densityWeights) {
                            const weight = state.densityWeights[i] || 0.5;
                            brightnessMultiplier = maxBrightnessMultiplier * weight * 1.5;
                        }
                        
                        // 亮度范围从 minBrightness 到 brightnessMultiplier
                        const brightness = (minBrightness + individualWave * (brightnessMultiplier - minBrightness)) * sparkle;
                        
                        const i3 = i * 3;
                        // 基于原始颜色应用亮度调整
                        colorArray[i3] = Math.min(1.0, state.baseColors[i3] * brightness);
                        colorArray[i3 + 1] = Math.min(1.0, state.baseColors[i3 + 1] * brightness);
                        colorArray[i3 + 2] = Math.min(1.0, state.baseColors[i3 + 2] * brightness);
                    }
                    colors.needsUpdate = true;
                }
            }
        } else if (state.heartbeat.mode === 'breath') {
            // 呼吸模式：平滑的正弦波呼吸效果
            state.heartbeat.phase += 0.012 * frequencyScale;
            if (state.heartbeat.phase > Math.PI * 2) {
                state.heartbeat.phase -= Math.PI * 2;
            }
            
            // 平滑的呼吸波形 (使用 cos 使得从亮开始)
            const breathWave = 0.5 + 0.5 * Math.cos(state.heartbeat.phase);
            
            // 最大亮度调整：100 + (n-100)*0.1
            const rawBrightness = state.heartbeat.maxBrightness || 200;
            const adjustedBrightness = 100 + (rawBrightness - 100) * 0.1;
            const maxStrength = adjustedBrightness / 100 * 1.5 * 1.5;
            state.bloomPass.strength = breathWave * maxStrength;
            
            const minOpacity = 0.25;
            const targetOpacity = minOpacity + breathWave * (1.0 - minOpacity);
            
            if (state.points && state.points.material) {
                state.points.material.opacity = targetOpacity;
                state.points.material.needsUpdate = true;
            }
        } else if (state.heartbeat.mode === 'diffuse') {
            // 扩散模式：从高密度区域向低密度区域扩散的波纹
            state.heartbeat.phase += 0.02 * frequencyScale;
            
            // 最大亮度调整：100 + (n-100)*0.1
            const rawBrightness = state.heartbeat.maxBrightness || 200;
            const adjustedBrightness = 100 + (rawBrightness - 100) * 0.1;
            const maxStrength = adjustedBrightness / 100 * 1.5 * 1.5;
            state.bloomPass.strength = maxStrength * 0.7;
            
            if (state.points && state.points.geometry) {
                const colors = state.points.geometry.attributes.color;
                const positions = state.points.geometry.attributes.position;
                if (colors && positions) {
                    const pointCount = state.pointCount || 0;
                    const time = state.heartbeat.phase;
                    
                    const maxPoints = state.MAX_POINTS || 100000;
                    if (!state.baseColors || state.baseColors.length < maxPoints * 3) {
                        state.baseColors = new Float32Array(maxPoints * 3);
                        state.baseColorsCount = 0;
                    }
                    if (state.baseColorsCount === 0 && pointCount > 0) {
                        const colorArray = colors.array;
                        for (let i = 0; i < pointCount * 3; i++) {
                            state.baseColors[i] = colorArray[i];
                        }
                        state.baseColorsCount = pointCount;
                    }
                    
                    // 计算概率密度用于扩散效果
                    const orbitalKey = state.currentOrbital || '1s';
                    const Hydrogen = window.Hydrogen;
                    let orbitalParams = null;
                    if (Hydrogen && Hydrogen.orbitalParamsFromKey) {
                        orbitalParams = Hydrogen.orbitalParamsFromKey(orbitalKey);
                    }
                    
                    // 初始化密度数组
                    if (!state.diffuseDensities || state.diffuseDensities.length < maxPoints) {
                        state.diffuseDensities = new Float32Array(maxPoints);
                        state.diffuseDensitiesComputed = false;
                    }
                    
                    const colorArray = colors.array;
                    const processCount = Math.min(pointCount, state.baseColorsCount);
                    
                    // 如果密度未计算，先计算一次
                    if (!state.diffuseDensitiesComputed && orbitalParams && Hydrogen.density3D_real) {
                        let maxDensity = 0;
                        for (let i = 0; i < processCount; i++) {
                            const i3 = i * 3;
                            const x = positions.array[i3];
                            const y = positions.array[i3 + 1];
                            const z = positions.array[i3 + 2];
                            const r = Math.sqrt(x * x + y * y + z * z);
                            const theta = r > 0 ? Math.acos(z / r) : 0;
                            const phi = Math.atan2(y, x);
                            
                            const density = Hydrogen.density3D_real(
                                orbitalParams.angKey,
                                orbitalParams.n,
                                orbitalParams.l,
                                r, theta, phi
                            );
                            state.diffuseDensities[i] = density;
                            if (density > maxDensity) maxDensity = density;
                        }
                        // 归一化密度到 0-1 范围
                        if (maxDensity > 0) {
                            for (let i = 0; i < processCount; i++) {
                                state.diffuseDensities[i] = state.diffuseDensities[i] / maxDensity;
                            }
                        }
                        state.diffuseDensitiesComputed = true;
                    }
                    
                    for (let i = 0; i < processCount; i++) {
                        const i3 = i * 3;
                        
                        // 使用密度值作为扩散相位基准
                        // 密度高的点先亮，然后波浪向密度低的区域扩散
                        const density = state.diffuseDensities[i] || 0;
                        // 用 log 压缩使效果更明显
                        const logDensity = Math.log(density * 99 + 1) / Math.log(100);
                        
                        // 波浪从高密度向低密度传播
                        const wavePhase = (1 - logDensity) * 6 - time * 3;
                        const wave = 0.5 + 0.5 * Math.sin(wavePhase);
                        const brightness = 0.3 + wave * 1.2;
                        
                        colorArray[i3] = Math.min(1.0, state.baseColors[i3] * brightness);
                        colorArray[i3 + 1] = Math.min(1.0, state.baseColors[i3 + 1] * brightness);
                        colorArray[i3 + 2] = Math.min(1.0, state.baseColors[i3 + 2] * brightness);
                    }
                    colors.needsUpdate = true;
                }
            }
        } else if (state.heartbeat.mode === 'wave') {
            // 波浪模式：沿等值面扩散 - 点亮等值面轮廓，从内向外扩散
            // 等值面是包含特定百分比点的表面，比如P轨道的哑铃形
            state.heartbeat.phase += 0.00375 * frequencyScale; // 原来0.015，现在缩小到1/4
            
            // 最大亮度调整：100 + (n-100)*0.1
            const rawBrightness = state.heartbeat.maxBrightness || 200;
            const adjustedBrightness = 100 + (rawBrightness - 100) * 0.1;
            const maxStrength = adjustedBrightness / 100 * 1.5 * 1.5;
            state.bloomPass.strength = maxStrength * 0.9;
            
            if (state.points && state.points.geometry) {
                const colors = state.points.geometry.attributes.color;
                const positions = state.points.geometry.attributes.position;
                if (colors && positions) {
                    const pointCount = state.pointCount || 0;
                    const time = state.heartbeat.phase;
                    
                    const maxPoints = state.MAX_POINTS || 100000;
                    if (!state.baseColors || state.baseColors.length < maxPoints * 3) {
                        state.baseColors = new Float32Array(maxPoints * 3);
                        state.baseColorsCount = 0;
                    }
                    if (state.baseColorsCount === 0 && pointCount > 0) {
                        const colorArray = colors.array;
                        for (let i = 0; i < pointCount * 3; i++) {
                            state.baseColors[i] = colorArray[i];
                        }
                        state.baseColorsCount = pointCount;
                    }
                    
                    const Hydrogen = window.Hydrogen;
                    const ui = window.ElectronCloud.ui;
                    
                    // 检查是否为比照模式（多轨道独立显示）
                    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
                    const orbitals = state.currentOrbitals || [state.currentOrbital || '1s'];
                    
                    // 缓存所有轨道的参数
                    if (!state.waveOrbitalParams) {
                        state.waveOrbitalParams = [];
                    }
                    
                    // 初始化数组
                    if (!state.waveRanks || state.waveRanks.length < maxPoints) {
                        state.waveRanks = new Float32Array(maxPoints); // 每个点的等值面排名(0-1)
                        state.waveDensities = new Float32Array(maxPoints);
                        state.waveRanksComputed = false;
                        state.waveRanksPointCount = 0;
                    }
                    
                    const colorArray = colors.array;
                    const processCount = Math.min(pointCount, state.baseColorsCount);
                    
                    // 当点数变化时，重新计算等值面排名
                    // 每隔 5000 点重新计算一次，避免频繁更新影响性能
                    const shouldRecompute = !state.waveRanksComputed || 
                        (processCount >= 5000 && Math.floor(processCount / 5000) > Math.floor(state.waveRanksPointCount / 5000));
                    
                    // 计算等值面排名：基于密度值排序，得到每个点属于哪个等值面
                    if (shouldRecompute && Hydrogen && Hydrogen.density3D_real && Hydrogen.orbitalParamsFromKey) {
                        // 更新轨道参数缓存
                        state.waveOrbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
                        
                        if (state.waveOrbitalParams.length > 0) {
                            const usePerPointOrbital = isCompareMode && state.pointOrbitalIndices && state.waveOrbitalParams.length > 1;
                            
                            // 第一步：计算每个点的密度
                            for (let i = 0; i < processCount; i++) {
                                const i3 = i * 3;
                                const x = positions.array[i3];
                                const y = positions.array[i3 + 1];
                                const z = positions.array[i3 + 2];
                                const r = Math.sqrt(x * x + y * y + z * z);
                                const theta = r > 0 ? Math.acos(z / r) : 0;
                                const phi = Math.atan2(y, x);
                                
                                // 根据模式选择轨道参数和密度计算方式
                                let density = 0;
                                if (usePerPointOrbital) {
                                    // 比照模式：使用每个点实际所属的轨道
                                    const orbitalIndex = state.pointOrbitalIndices[i] || 0;
                                    const orbitalParams = state.waveOrbitalParams[orbitalIndex] || state.waveOrbitalParams[0];
                                    density = Hydrogen.density3D_real(
                                        orbitalParams.angKey,
                                        orbitalParams.n,
                                        orbitalParams.l,
                                        r, theta, phi
                                    );
                                } else {
                                    // 【关键修复】多选叠加模式：使用所有轨道的混合密度
                                    // 等值面应该基于叠加后的概率密度分布
                                    for (const orbitalParams of state.waveOrbitalParams) {
                                        density += Hydrogen.density3D_real(
                                            orbitalParams.angKey,
                                            orbitalParams.n,
                                            orbitalParams.l,
                                            r, theta, phi
                                        );
                                    }
                                    // 取平均（与采样和理论曲线计算一致）
                                    density /= state.waveOrbitalParams.length;
                                }
                                
                                state.waveDensities[i] = density;
                            }
                            
                            if (usePerPointOrbital) {
                                // 比照模式：对每个轨道的点分别排序
                                // 这样每个轨道的等值面是独立计算的
                                const orbitalPointGroups = {};
                                
                                // 按轨道分组
                                for (let i = 0; i < processCount; i++) {
                                    const orbitalIndex = state.pointOrbitalIndices[i] || 0;
                                    if (!orbitalPointGroups[orbitalIndex]) {
                                        orbitalPointGroups[orbitalIndex] = [];
                                    }
                                    orbitalPointGroups[orbitalIndex].push(i);
                                }
                                
                                // 对每个轨道组内部按密度排序
                                for (const groupKey of Object.keys(orbitalPointGroups)) {
                                    const groupIndices = orbitalPointGroups[groupKey];
                                    groupIndices.sort((a, b) => state.waveDensities[b] - state.waveDensities[a]);
                                    
                                    // 计算组内排名
                                    const groupSize = groupIndices.length;
                                    for (let rank = 0; rank < groupSize; rank++) {
                                        const pointIndex = groupIndices[rank];
                                        state.waveRanks[pointIndex] = rank / (groupSize - 1 || 1);
                                    }
                                }
                            } else {
                                // 正常模式：所有点一起排序
                                const indices = new Array(processCount);
                                for (let i = 0; i < processCount; i++) indices[i] = i;
                                indices.sort((a, b) => state.waveDensities[b] - state.waveDensities[a]);
                                
                                // 计算每个点的等值面排名
                                // 排名 0 = 最高密度（最内层等值面，约0%）
                                // 排名 1 = 最低密度（最外层等值面，约99%）
                                for (let rank = 0; rank < processCount; rank++) {
                                    const pointIndex = indices[rank];
                                    state.waveRanks[pointIndex] = rank / (processCount - 1 || 1);
                                }
                            }
                            
                            state.waveRanksComputed = true;
                            state.waveRanksPointCount = processCount;
                        }
                    }
                    
                    // 波浪扩散：单向从核心向外扩散
                    // 使用 smoothstep 曲线让周期衔接自然
                    const cycleProgress = time % 1.0; // 周期为1
                    
                    // smoothstep: 3t² - 2t³
                    // 开始时慢（等值面从-0.05缓慢启动，此时全黑）
                    // 中间快（快速穿过可见区域 0~1）
                    // 结束时慢（等值面缓慢到达1.05，此时全黑）
                    const t = cycleProgress;
                    const smoothT = t * t * (3 - 2 * t);
                    
                    // 等值面位置：从 -0.05 到 1.05
                    // 起点和终点都超出范围，确保开始和结束时全黑
                    const pulsePosition = -0.05 + smoothT * 1.1;
                    
                    // 等值面薄层发光
                    const pulseWidth = 0.025;
                    
                    for (let i = 0; i < processCount; i++) {
                        const i3 = i * 3;
                        
                        // 点的等值面排名 (0=最高密度/内层, 1=最低密度/外层)
                        const rank = state.waveRanks[i] || 0;
                        
                        // 计算距离等值面的距离
                        const distToPulse = Math.abs(rank - pulsePosition);
                        
                        // 薄层发光：只有非常接近等值面的点才亮，其他全黑
                        let brightness;
                        if (distToPulse < pulseWidth) {
                            // 核心区域：最亮
                            brightness = 2.5;
                        } else if (distToPulse < pulseWidth * 2) {
                            // 边缘区域：平滑衰减
                            const fadeT = (distToPulse - pulseWidth) / pulseWidth;
                            brightness = 2.5 * (1 - fadeT * fadeT);
                        } else {
                            // 远离等值面：全黑
                            brightness = 0;
                        }
                        
                        colorArray[i3] = Math.min(1.0, state.baseColors[i3] * brightness);
                        colorArray[i3 + 1] = Math.min(1.0, state.baseColors[i3 + 1] * brightness);
                        colorArray[i3 + 2] = Math.min(1.0, state.baseColors[i3 + 2] * brightness);
                    }
                    colors.needsUpdate = true;
                }
            }
        }
        state.bloomEnabled = true;
    }
    
    // 独立权重模式：即使闪烁模式关闭，也可以根据密度权重调整颜色亮度
    // 这使得静态辉光也能按权重显示
    // 辉光亮度 = 基础亮度 × 权重系数(不开启则为1) × 闪烁系数(不开启则为1)
    if (!state.heartbeat.enabled && state.points && state.points.geometry) {
        // 非闪烁模式下，始终启用辉光，只是强度不同
        state.bloomEnabled = true;
        
        const colors = state.points.geometry.attributes.color;
        const positions = state.points.geometry.attributes.position;
        if (colors && positions) {
            const pointCount = state.pointCount || 0;
            const maxPoints = state.MAX_POINTS || 100000;
            
            // 确保有基础颜色数组
            if (!state.baseColors || state.baseColors.length < maxPoints * 3) {
                state.baseColors = new Float32Array(maxPoints * 3);
                state.baseColorsCount = 0;
            }
            
            // 保存原始颜色（如果还没有）
            if (state.baseColorsCount === 0 && pointCount > 0) {
                const colorArray = colors.array;
                for (let i = 0; i < pointCount * 3; i++) {
                    state.baseColors[i] = colorArray[i];
                }
                state.baseColorsCount = pointCount;
            }
            
            const colorArray = colors.array;
            const processCount = Math.min(pointCount, state.baseColorsCount);
            
            if (state.weightMode) {
                // 权重模式开启：启用bloom，设置较高的辉光强度
                state.bloomEnabled = true;
                if (state.bloomPass) {
                    state.bloomPass.strength = 1.5; // 权重模式辉光强度
                }
                
                // 初始化/更新密度权重
                if (!state.densityWeights || state.densityWeights.length < maxPoints) {
                    state.densityWeights = new Float32Array(maxPoints);
                    for (let i = 0; i < pointCount; i++) {
                        const i3 = i * 3;
                        const x = positions.array[i3];
                        const y = positions.array[i3 + 1];
                        const z = positions.array[i3 + 2];
                        const r = Math.sqrt(x * x + y * y + z * z);
                        state.densityWeights[i] = Math.exp(-r / 10) * 0.7 + 0.3;
                    }
                }
                
                // 根据权重调整颜色亮度
                for (let i = 0; i < processCount; i++) {
                    const weight = state.densityWeights[i] || 0.5;
                    const brightness = 0.5 + weight * 1.0; // 0.5 - 1.5 范围
                    const i3 = i * 3;
                    colorArray[i3] = Math.min(1.0, state.baseColors[i3] * brightness);
                    colorArray[i3 + 1] = Math.min(1.0, state.baseColors[i3 + 1] * brightness);
                    colorArray[i3 + 2] = Math.min(1.0, state.baseColors[i3 + 2] * brightness);
                }
                colors.needsUpdate = true;
            } else {
                // 权重模式关闭：恢复原始颜色（权重系数 = 1）
                // 【修复】使用亮度滑动条设置的值，而不是固定像素
                if (state.bloomPass && state.bloomStrength !== undefined) {
                    state.bloomPass.strength = state.bloomStrength;
                }
                
                // 只在颜色被修改过的情况下才恢复
                if (state.weightModeWasEnabled) {
                    for (let i = 0; i < processCount; i++) {
                        const i3 = i * 3;
                        colorArray[i3] = state.baseColors[i3];
                        colorArray[i3 + 1] = state.baseColors[i3 + 1];
                        colorArray[i3 + 2] = state.baseColors[i3 + 2];
                    }
                    colors.needsUpdate = true;
                    state.weightModeWasEnabled = false;
                }
            }
            
            // 记录权重模式状态，用于下一帧检测关闭
            if (state.weightMode) {
                state.weightModeWasEnabled = true;
            }
        }
    }

    state.controls.update();
    
    // 根据辉光状态选择渲染方式
    if (state.bloomEnabled && state.composer) {
        // 先隐藏坐标轴（不参与bloom）
        const axesWasVisible = state.customAxes ? state.customAxes.visible : false;
        const axisHelperWasVisible = state.autoRotate.axisHelper ? state.autoRotate.axisHelper.visible : false;
        
        if (state.customAxes) state.customAxes.visible = false;
        if (state.autoRotate.axisHelper) state.autoRotate.axisHelper.visible = false;
        
        // 隐藏点云，只做bloom渲染
        const pointsWasVisible = state.points ? state.points.visible : false;
        const angularWasVisible = state.angularOverlay ? state.angularOverlay.visible : false;
        
        // 渲染bloom效果（只有点云参与）
        state.composer.render();
        
        // 恢复坐标轴可见性
        if (state.customAxes) state.customAxes.visible = axesWasVisible;
        if (state.autoRotate.axisHelper) state.autoRotate.axisHelper.visible = axisHelperWasVisible;
        
        // 单独渲染坐标轴（无bloom，叠加在上面）
        if (axesWasVisible || axisHelperWasVisible) {
            // 临时隐藏点云，只渲染坐标轴
            if (state.points) state.points.visible = false;
            if (state.angularOverlay) state.angularOverlay.visible = false;
            
            state.renderer.autoClear = false;
            state.renderer.clearDepth();
            state.renderer.render(state.scene, state.camera);
            state.renderer.autoClear = true;
            
            // 恢复点云可见性
            if (state.points) state.points.visible = pointsWasVisible;
            if (state.angularOverlay) state.angularOverlay.visible = angularWasVisible;
        }
    } else {
        state.renderer.render(state.scene, state.camera);
    }
};

// 采样完成时的处理
window.ElectronCloud.Scene.onSamplingCompleted = function() {
    const state = window.ElectronCloud.state;
    
    state.samplingCompleted = true;
    state.renderingCompleted = true;
    
    // 更新坐标轴滑动条状态（采样完成后可修改）
    window.ElectronCloud.UI.updateAxesSizeRangeState();
    
    // 更新自动旋转按钮状态（采样完成后可启用）
    if (window.ElectronCloud.UI.updateAutoRotateButtonState) {
        window.ElectronCloud.UI.updateAutoRotateButtonState();
    }
    
    console.log('采样完成，准备显示数据面板');
    
    // 初始化所有轨道为可见状态
    for (const orbitalKey of state.currentOrbitals) {
        state.orbitalVisibility[orbitalKey] = true;
    }
    
    // 如果在比照模式下，更新标签提示
    const compareToggle = window.ElectronCloud.ui.compareToggle;
    if (compareToggle && compareToggle.checked) {
        // 渲染完成后更新标签提示和选择计数
        window.ElectronCloud.UI.updateSelectionCount();
        
        // 添加渲染完成状态的CSS类，用于hover效果
        const controlPanel = document.getElementById('control-panel');
        if (controlPanel) {
            controlPanel.classList.add('rendering-completed');
        }
    }
    
    // 最后一次更新后台数据
    window.ElectronCloud.Orbital.updateBackgroundChartData();
    
    // 采样完成后确保数据面板可见（不改变用户的展开/收起状态）
    const dataPanel = document.getElementById('data-panel');
    if (dataPanel) {
        dataPanel.classList.remove('collapsed');
    }
    
    // 采样完成后立即显示图表
    window.ElectronCloud.Orbital.drawProbabilityChart(true);
    
    // 停止绘制
    state.isDrawing = false;
    
    // 更新3D角向开关状态（采样完成后应该可用）
    window.ElectronCloud.UI.updateAngular3DToggleState();

    // 更新全屏按钮状态
    window.ElectronCloud.UI.updateFullscreenBtnState();

    // 启用手势按钮
    window.ElectronCloud.UI.enableGestureButton();
};

// 清除场景中的点云
window.ElectronCloud.Scene.clearPoints = function() {
    const state = window.ElectronCloud.state;
    
    if (state.points) {
        state.scene.remove(state.points);
        state.points.geometry.dispose();
        state.points.material.dispose();
        state.points = null;
    }
};

// 清除角向分布叠加
window.ElectronCloud.Scene.clearAngularOverlay = function() {
    const state = window.ElectronCloud.state;
    
    if (state.angularOverlay) {
        state.scene.remove(state.angularOverlay);
        state.angularOverlay.traverse((child) => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
        state.angularOverlay = null;
    }
};

// 【根本修复】统一重置所有场景对象的旋转状态
// 这是解决"坐标歪斜"问题的关键函数
// 问题根本原因：points、customAxes、angularOverlay 三个对象的旋转状态可能不同步
// 当新建点云时，旧的坐标轴可能保持之前的旋转角度，导致不对齐
window.ElectronCloud.Scene.resetAllSceneObjectsRotation = function() {
    const state = window.ElectronCloud.state;
    
    console.log('【统一旋转重置】重置所有场景对象的旋转状态');
    
    // 使用UI模块的公共函数重置旋转（避免代码重复）
    if (window.ElectronCloud.UI && window.ElectronCloud.UI.resetSceneObjectsRotation) {
        window.ElectronCloud.UI.resetSceneObjectsRotation();
    } else {
        // 降级方案：直接重置（防止初始化顺序问题）
        if (state.points) {
            state.points.rotation.set(0, 0, 0);
            state.points.updateMatrix();
            state.points.updateMatrixWorld(true);
        }
        if (state.angularOverlay) {
            state.angularOverlay.rotation.set(0, 0, 0);
            state.angularOverlay.updateMatrix();
            state.angularOverlay.updateMatrixWorld(true);
        }
        if (state.customAxes) {
            state.customAxes.rotation.set(0, 0, 0);
            state.customAxes.updateMatrix();
            state.customAxes.updateMatrixWorld(true);
        }
    }
    
    // 重置旋转轴辅助线（如果存在）
    if (state.autoRotate && state.autoRotate.axisHelper) {
        state.scene.remove(state.autoRotate.axisHelper);
        state.autoRotate.axisHelper.geometry.dispose();
        state.autoRotate.axisHelper.material.dispose();
        state.autoRotate.axisHelper = null;
    }
    
    // 重置自动旋转累计角度
    if (state.autoRotate) {
        state.autoRotate.totalAngle = 0;
    }
    
    // 如果有锁定视角，清除锁定状态
    if (state.lockedAxis) {
        state.lockedAxis = null;
        // 恢复OrbitControls（但保持 enableRotate = false，使用自定义轨迹球）
        if (state.controls) {
            state.controls.enabled = true;
            // 注意：不恢复 enableRotate，因为使用自定义轨迹球旋转
            state.controls.enableZoom = true;
            // 检查 centerLock 状态，决定是否启用平移
            const centerLock = document.getElementById('center-lock');
            state.controls.enablePan = !(centerLock && centerLock.checked);
        }
        // 恢复相机默认up向量
        if (state.camera) {
            state.camera.up.set(0, 1, 0);
        }
        // 恢复所有坐标轴的可见性
        if (window.ElectronCloud.Scene.resetAxesVisibility) {
            window.ElectronCloud.Scene.resetAxesVisibility();
        }
        // 清除锁定轴旋转处理器
        if (window.ElectronCloud.UI.clearLockedAxisRotation) {
            window.ElectronCloud.UI.clearLockedAxisRotation();
        }
    }
    
    console.log('【统一旋转重置】完成');
};

// 生成圆形点的 sprite 纹理（带缓存）
// 【优化】使用更锐利的边缘，减少摩尔纹干涉
let _circleSpriteCached = null;
window.ElectronCloud.Scene.generateCircleSprite = function() {
    // 如果已有缓存，直接返回
    if (_circleSpriteCached) return _circleSpriteCached;
    
    const size = 64;
    const canvas = document.createElement('canvas');
    canvas.width = size;
    canvas.height = size;
    const ctx = canvas.getContext('2d');

    // 背景透明
    ctx.clearRect(0, 0, size, size);

    // 【柔和边缘】使用平滑的高斯型边缘过渡，减少视觉干扰
    const grad = ctx.createRadialGradient(size/2, size/2, 0, size/2, size/2, size/2);
    grad.addColorStop(0, 'rgba(255,255,255,1)');
    grad.addColorStop(0.2, 'rgba(255,255,255,0.9)');
    grad.addColorStop(0.4, 'rgba(255,255,255,0.6)');
    grad.addColorStop(0.6, 'rgba(255,255,255,0.3)');
    grad.addColorStop(0.8, 'rgba(255,255,255,0.1)');
    grad.addColorStop(1, 'rgba(255,255,255,0)');
    ctx.fillStyle = grad;
    ctx.beginPath();
    ctx.arc(size/2, size/2, size/2.2, 0, Math.PI * 2);
    ctx.fill();

    const texture = new THREE.CanvasTexture(canvas);
    texture.minFilter = THREE.LinearFilter;
    
    // 缓存纹理
    _circleSpriteCached = texture;
    return texture;
};

// 显示旋转轴辅助线（虚线，1秒后消失）
window.ElectronCloud.Scene.showRotationAxisHelper = function(axis) {
    const state = window.ElectronCloud.state;
    
    // 移除旧的辅助线
    if (state.autoRotate.axisHelper) {
        state.scene.remove(state.autoRotate.axisHelper);
        state.autoRotate.axisHelper.geometry.dispose();
        state.autoRotate.axisHelper.material.dispose();
        state.autoRotate.axisHelper = null;
    }
    
    // 归一化轴向量
    const normalizedAxis = axis.clone().normalize();
    if (normalizedAxis.length() === 0) return;
    
    // 计算线段长度（基于当前视图范围）
    const lineLength = Math.max(50, state.farthestDistance * 1.5);
    
    // 创建虚线几何体
    const points = [
        normalizedAxis.clone().multiplyScalar(-lineLength),
        normalizedAxis.clone().multiplyScalar(lineLength)
    ];
    const geometry = new THREE.BufferGeometry().setFromPoints(points);
    
    // 虚线材质（红色）
    const material = new THREE.LineDashedMaterial({
        color: 0xff3333,
        dashSize: 2,
        gapSize: 1,
        linewidth: 2
    });
    
    const line = new THREE.Line(geometry, material);
    line.computeLineDistances(); // 虚线需要这个
    
    // 放在layer 1，不受bloom影响
    line.layers.set(1);
    
    state.autoRotate.axisHelper = line;
    state.scene.add(line);
    
    // 1秒后自动消失
    setTimeout(() => {
        if (state.autoRotate.axisHelper === line) {
            state.scene.remove(line);
            line.geometry.dispose();
            line.material.dispose();
            state.autoRotate.axisHelper = null;
        }
    }, 1000);
};

// ========================================
// 自定义轨迹球鼠标控制（完全自由旋转，无极点限制，带惯性）
// ========================================
window.ElectronCloud.Scene.initTrackballMouseControl = function() {
    const state = window.ElectronCloud.state;
    const canvas = state.renderer.domElement;
    
    // 配置参数 - 优化手感，更接近原 OrbitControls
    const CONFIG = {
        rotationSpeed: 0.005,      // 旋转灵敏度（稍微降低，更精细）
        friction: 0.92,            // 惯性摩擦系数 (0-1，降低使惯性衰减更快)
        minVelocity: 0.0001,       // 停止惯性的最小速度（提高阈值，更快停止）
        velocitySamples: 3,        // 用于计算释放速度的样本数（减少延迟感）
        inertiaMultiplier: 0.8,    // 惯性放大倍数（降低，惯性更自然）
    };
    
    // 状态变量
    let isDragging = false;
    let previousMouseX = 0;
    let previousMouseY = 0;
    let activePointerId = null;
    
    // 惯性相关
    let velocityX = 0;
    let velocityY = 0;
    let velocityHistory = [];      // 记录最近几帧的速度
    let lastMoveTime = 0;
    let inertiaAnimationId = null;
    
    // 轨迹球旋转实现
    function applyTrackballRotation(deltaX, deltaY) {
        // 如果视角被锁定或手势控制正在运行，不处理
        if (state.lockedAxis) return;
        if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.isRunning && window.ElectronCloud.Gesture.isRunning()) return;
        
        const camera = state.camera;
        const controls = state.controls;
        const target = controls.target;

        // 计算相机相对于目标点的偏移
        const offset = new THREE.Vector3().subVectors(camera.position, target);
        
        // 获取当前相机的坐标系
        const cameraDirection = new THREE.Vector3();
        camera.getWorldDirection(cameraDirection);
        
        // 计算相机的右轴和上轴（基于当前相机姿态）
        const cameraRight = new THREE.Vector3();
        cameraRight.crossVectors(cameraDirection, camera.up).normalize();
        
        // 如果 cameraRight 接近零（相机 up 与观察方向平行），使用备用方案
        if (cameraRight.lengthSq() < 0.0001) {
            if (Math.abs(cameraDirection.x) < 0.9) {
                cameraRight.crossVectors(new THREE.Vector3(1, 0, 0), cameraDirection).normalize();
            } else {
                cameraRight.crossVectors(new THREE.Vector3(0, 1, 0), cameraDirection).normalize();
            }
        }
        
        // 相机的实际上轴（与观察方向正交）
        const cameraUp = new THREE.Vector3();
        cameraUp.crossVectors(cameraRight, cameraDirection).normalize();
        
        // 创建旋转四元数
        const qHorizontal = new THREE.Quaternion().setFromAxisAngle(cameraUp, -deltaX);
        const qVertical = new THREE.Quaternion().setFromAxisAngle(cameraRight, -deltaY);
        
        // 合并旋转
        const quaternion = new THREE.Quaternion();
        quaternion.multiplyQuaternions(qHorizontal, qVertical);
        
        // 应用旋转到偏移向量
        offset.applyQuaternion(quaternion);
        
        // 更新相机位置
        camera.position.copy(target).add(offset);
        
        // 同时旋转 camera.up 向量，实现完全自由旋转
        camera.up.applyQuaternion(quaternion).normalize();
        
        // 防止 up 向量数值漂移（正交化校正）
        const newDirection = new THREE.Vector3().subVectors(target, camera.position).normalize();
        const dotProduct = camera.up.dot(newDirection);
        if (Math.abs(dotProduct) > 0.0001) {
            camera.up.sub(newDirection.clone().multiplyScalar(dotProduct)).normalize();
        }
        
        // 让相机看向目标点
        camera.lookAt(target);
        
        // 更新控制器目标
        controls.target.copy(target);
    }
    
    // 惯性动画循环
    function inertiaLoop() {
        // 如果正在拖动或视角被锁定，停止惯性
        if (isDragging || state.lockedAxis) {
            velocityX = 0;
            velocityY = 0;
            return;
        }
        
        // 如果手势控制正在运行，停止惯性
        if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.isRunning && window.ElectronCloud.Gesture.isRunning()) {
            velocityX = 0;
            velocityY = 0;
            return;
        }
        
        // 检查速度是否足够小可以停止
        const speed = Math.sqrt(velocityX * velocityX + velocityY * velocityY);
        if (speed < CONFIG.minVelocity) {
            velocityX = 0;
            velocityY = 0;
            return;
        }
        
        // 应用旋转
        applyTrackballRotation(velocityX, velocityY);
        
        // 应用摩擦力衰减
        velocityX *= CONFIG.friction;
        velocityY *= CONFIG.friction;
        
        // 继续动画
        inertiaAnimationId = requestAnimationFrame(inertiaLoop);
    }
    
    // 指针按下
    function onPointerDown(e) {
        // 只响应左键（鼠标）
        if (e.pointerType === 'mouse' && e.button !== 0) return;
        
        // 如果视角被锁定，不处理
        if (state.lockedAxis) return;
        
        // 如果手势控制正在运行，不处理
        if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.isRunning && window.ElectronCloud.Gesture.isRunning()) return;
        
        // 只处理鼠标左键的旋转
        if (e.pointerType === 'mouse' && e.button === 0) {
            isDragging = true;
            activePointerId = e.pointerId;
            previousMouseX = e.clientX;
            previousMouseY = e.clientY;
            lastMoveTime = performance.now();
            
            // 停止惯性
            velocityX = 0;
            velocityY = 0;
            velocityHistory = [];
            if (inertiaAnimationId) {
                cancelAnimationFrame(inertiaAnimationId);
                inertiaAnimationId = null;
            }
            
            // 捕获指针
            canvas.setPointerCapture(e.pointerId);
        }
    }
    
    // 指针移动
    function onPointerMove(e) {
        if (!isDragging || e.pointerId !== activePointerId) return;
        
        const currentTime = performance.now();
        const dt = currentTime - lastMoveTime;
        
        const deltaX = (e.clientX - previousMouseX) * CONFIG.rotationSpeed;
        const deltaY = (e.clientY - previousMouseY) * CONFIG.rotationSpeed;
        
        previousMouseX = e.clientX;
        previousMouseY = e.clientY;
        lastMoveTime = currentTime;
        
        // 记录速度历史（用于计算释放时的惯性）
        if (dt > 0) {
            velocityHistory.push({ vx: deltaX, vy: deltaY, dt: dt });
            // 只保留最近几个样本
            if (velocityHistory.length > CONFIG.velocitySamples) {
                velocityHistory.shift();
            }
        }
        
        applyTrackballRotation(deltaX, deltaY);
    }
    
    // 指针释放
    function onPointerUp(e) {
        if (e.pointerId === activePointerId) {
            isDragging = false;
            activePointerId = null;
            
            try {
                canvas.releasePointerCapture(e.pointerId);
            } catch (err) {
                // 忽略
            }
            
            // 计算释放时的速度（基于最近几帧的平均值）
            if (velocityHistory.length > 0) {
                let totalVx = 0;
                let totalVy = 0;
                let totalWeight = 0;
                
                // 加权平均，最近的帧权重更大
                velocityHistory.forEach((sample, index) => {
                    const weight = index + 1;  // 越新的样本权重越大
                    totalVx += sample.vx * weight;
                    totalVy += sample.vy * weight;
                    totalWeight += weight;
                });
                
                if (totalWeight > 0) {
                    velocityX = (totalVx / totalWeight) * CONFIG.inertiaMultiplier;
                    velocityY = (totalVy / totalWeight) * CONFIG.inertiaMultiplier;
                    
                    // 启动惯性动画
                    if (Math.abs(velocityX) > CONFIG.minVelocity || Math.abs(velocityY) > CONFIG.minVelocity) {
                        inertiaAnimationId = requestAnimationFrame(inertiaLoop);
                    }
                }
            }
            
            velocityHistory = [];
        }
    }
    
    // 绑定事件 - 使用捕获阶段
    canvas.addEventListener('pointerdown', onPointerDown, true);
    canvas.addEventListener('pointermove', onPointerMove, true);
    canvas.addEventListener('pointerup', onPointerUp, true);
    canvas.addEventListener('pointercancel', onPointerUp, true);
    
    // 保存引用
    state.trackballMouseControl = {
        enabled: true,
        handlers: {
            pointerdown: onPointerDown,
            pointermove: onPointerMove,
            pointerup: onPointerUp
        },
        // 暴露停止惯性的方法
        stopInertia: function() {
            velocityX = 0;
            velocityY = 0;
            if (inertiaAnimationId) {
                cancelAnimationFrame(inertiaAnimationId);
                inertiaAnimationId = null;
            }
        }
    };
    
    console.log('[Scene] 轨迹球鼠标控制已初始化（完全自由旋转 + 惯性）');
};

// 暂时禁用轨迹球鼠标控制
window.ElectronCloud.Scene.disableTrackballMouseControl = function() {
    const state = window.ElectronCloud.state;
    if (state.trackballMouseControl) {
        state.trackballMouseControl.enabled = false;
    }
};

// 恢复轨迹球鼠标控制
window.ElectronCloud.Scene.enableTrackballMouseControl = function() {
    const state = window.ElectronCloud.state;
    if (state.trackballMouseControl) {
        state.trackballMouseControl.enabled = true;
    }
};
