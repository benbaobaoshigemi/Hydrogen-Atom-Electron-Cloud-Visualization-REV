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

    // 渲染器
    const container = document.getElementById('container');
    state.renderer = new THREE.WebGLRenderer({ antialias: true });
    state.renderer.setSize(window.innerWidth, window.innerHeight);
    state.renderer.setPixelRatio(window.devicePixelRatio);
    container.appendChild(state.renderer.domElement);

    // 控制器
    state.controls = new THREE.OrbitControls(state.camera, state.renderer.domElement);
    state.controls.enableDamping = true;
    state.controls.dampingFactor = 0.05;

    // 自定义坐标轴（默认不显示）
    state.customAxes = window.ElectronCloud.Scene.createCustomAxes(constants.AXES_BASE_SIZE);
    state.customAxes.visible = false;
    state.scene.add(state.customAxes);

    // 监听窗口大小变化
    window.addEventListener('resize', window.ElectronCloud.Scene.onWindowResize, false);
    
    console.log('Three.js 场景初始化完成');
};

// 创建自定义坐标轴
window.ElectronCloud.Scene.createCustomAxes = function(size) {
    const axes = new THREE.Group();
    const material = new THREE.LineBasicMaterial({ color: 0xffffff });

    // 轴线
    const pointsX = [new THREE.Vector3(-size, 0, 0), new THREE.Vector3(size, 0, 0)];
    const pointsY = [new THREE.Vector3(0, -size, 0), new THREE.Vector3(0, size, 0)];
    const pointsZ = [new THREE.Vector3(0, 0, -size), new THREE.Vector3(0, 0, size)];
    
    axes.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsX), material));
    axes.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsY), material));
    axes.add(new THREE.Line(new THREE.BufferGeometry().setFromPoints(pointsZ), material));

    // 箭头
    const arrowGeom = new THREE.ConeGeometry(0.5, 2, 8);
    const arrowMat = new THREE.MeshBasicMaterial({ color: 0xffffff });

    const arrowX = new THREE.Mesh(arrowGeom, arrowMat);
    arrowX.position.set(size, 0, 0);
    arrowX.rotation.z = -Math.PI / 2;
    axes.add(arrowX);

    const arrowY = new THREE.Mesh(arrowGeom, arrowMat);
    arrowY.position.set(0, size, 0);
    axes.add(arrowY);

    const arrowZ = new THREE.Mesh(arrowGeom, arrowMat);
    arrowZ.position.set(0, 0, size);
    arrowZ.rotation.x = Math.PI / 2;
    axes.add(arrowZ);
    
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
            axes.add(textMeshX);

            const textShapesY = font.generateShapes('Y', 1.5);
            const textGeomY = new THREE.ShapeGeometry(textShapesY);
            const textMeshY = new THREE.Mesh(textGeomY, textMat);
            textMeshY.position.set(0, size + 2, 0);
            axes.add(textMeshY);

            const textShapesZ = font.generateShapes('Z', 1.5);
            const textGeomZ = new THREE.ShapeGeometry(textShapesZ);
            const textMeshZ = new THREE.Mesh(textGeomZ, textMat);
            textMeshZ.position.set(0, 0, size + 2);
            axes.add(textMeshZ);
        },
        undefined, // onProgress
        function (error) {
            console.warn('坐标轴字体加载失败，标签将不显示:', error);
        }
    );

    return axes;
};

// 处理窗口大小变化
window.ElectronCloud.Scene.onWindowResize = function() {
    const state = window.ElectronCloud.state;
    
    state.camera.aspect = window.innerWidth / window.innerHeight;
    state.camera.updateProjectionMatrix();
    state.renderer.setSize(window.innerWidth, window.innerHeight);
};

// 动画循环
window.ElectronCloud.Scene.animate = function() {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    
    state.animationFrameId = requestAnimationFrame(window.ElectronCloud.Scene.animate);

    // 如果正在绘制，则更新点
    if (state.isDrawing) {
        if (state.pointCount < state.MAX_POINTS) {
            window.ElectronCloud.Sampling.updatePoints();
            
            // 图表刷新：每秒刷新一次（1000ms），跳过第一秒
            const now = performance.now();
            if (!state.lastChartUpdateTime) state.lastChartUpdateTime = 0;
            
            // 检查是否已过1秒延迟期（只跳过第一次刷新）
            const isDelayPassed = state.samplingStartTime && (now - state.samplingStartTime > 1000);
            
            if (isDelayPassed && state.radialSamples.length > 0 && (now - state.lastChartUpdateTime) >= 1000) {
                state.lastChartUpdateTime = now;
                window.ElectronCloud.Orbital.updateBackgroundChartData();
                // 如果数据面板可见，实时刷新图表
                const dataPanel = document.getElementById('data-panel');
                if (dataPanel && !dataPanel.classList.contains('collapsed')) {
                    window.ElectronCloud.Orbital.drawProbabilityChart(false);
                }
            }
        }
        
        // 检查采样完成
        if (state.pointCount >= state.MAX_POINTS && !state.samplingCompleted) {
            window.ElectronCloud.Scene.onSamplingCompleted();
        }
    }

    // 确保坐标系状态和大小正确：实时更新坐标系
    if (state.customAxes && state.axesScaleFactor > 0) {
        if (state.farthestDistance === 0) {
            state.customAxes.visible = false;
        } else {
            // 当有了实际轨道半径且用户设置了比例系数时，显示坐标系
            state.customAxes.visible = true;
            
            // 实时更新坐标系大小：根据当前的farthestDistance和比例系数
            const orbitalRadius = Math.max(constants.AXES_BASE_SIZE, state.farthestDistance);
            const targetSize = orbitalRadius * state.axesScaleFactor;
            const scale = targetSize / constants.AXES_BASE_SIZE;
            
            // 只有当缩放比例发生变化时才更新（避免不必要的更新）
            const currentScale = state.customAxes.scale.x;
            if (Math.abs(currentScale - scale) > 0.001) {
                state.customAxes.scale.set(scale, scale, scale);
            }
        }
    } else if (state.customAxes) {
        // 比例系数为0时隐藏坐标系
        state.customAxes.visible = false;
    }

    state.controls.update();
    state.renderer.render(state.scene, state.camera);
};

// 采样完成时的处理
window.ElectronCloud.Scene.onSamplingCompleted = function() {
    const state = window.ElectronCloud.state;
    
    state.samplingCompleted = true;
    state.renderingCompleted = true;
    
    // 更新坐标轴滑动条状态（采样完成后可修改）
    window.ElectronCloud.UI.updateAxesSizeRangeState();
    
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

// 生成圆形点的 sprite 纹理（带缓存）
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

    // 白色圆形
    const grad = ctx.createRadialGradient(size/2, size/2, 0, size/2, size/2, size/2);
    grad.addColorStop(0, 'rgba(255,255,255,1)');
    grad.addColorStop(0.6, 'rgba(255,255,255,0.9)');
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
