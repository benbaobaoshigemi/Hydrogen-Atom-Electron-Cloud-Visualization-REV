// 核心变量和常量定义
window.ElectronCloud = window.ElectronCloud || {};

// 全局变量
window.ElectronCloud.state = {
    // Three.js 核心对象
    scene: null,
    camera: null,
    renderer: null,
    controls: null,
    
    // 渲染对象
    points: null,
    customAxes: null,
    angularOverlay: null,
    
    // 轨道状态
    currentOrbital: '1s',
    currentOrbitals: ['1s'],
    
    // 采样状态
    animationFrameId: null,
    pointCount: 0,
    MAX_POINTS: 50000,
    pointsPerFrame: 200,
    isDrawing: false,
    farthestDistance: 0,
    samplingBoundary: 40,
    samplingCompleted: false,
    roundRobinIndex: 0, // 【新增】多选模式轮流采样索引
    
    // 采样数据
    radialSamples: [],
    angularSamples: [],
    
    // 轨道显示状态管理
    orbitalVisibility: {},
    renderingCompleted: false,
    orbitalPointsMap: {},
    originalPositions: null,
    orbitalSamplesMap: {},
    pointOrbitalIndices: null, // 每个点的轨道索引（Int16Array）
    
    // 3D 相关
    angularUpdated: false,
    axesScaleFactor: 0, // 坐标轴比例系数 (0-1.05)
    
    // 图表相关
    probabilityChart: null,
    isLiveChart: true,
    backgroundChartData: {
        radial: null,
        angular: null
    }
};

// 常量
window.ElectronCloud.constants = {
    AXES_BASE_SIZE: 30,
    
    // 比照模式颜色方案
    compareColors: [
        { name: 'red', value: [1, 0.2, 0.2] },
        { name: 'green', value: [0.2, 1, 0.2] },
        { name: 'blue', value: [0.2, 0.2, 1] }
    ],
    
    // 轨道键值到显示名称的映射（统一定义，避免重复）
    orbitalDisplayNames: {
        '1s': '1s',
        '2s': '2s', '2px': '2px', '2py': '2py', '2pz': '2pz',
        '3s': '3s', '3px': '3px', '3py': '3py', '3pz': '3pz',
        '3d_xy': '3d(xy)', '3d_xz': '3d(xz)', '3d_yz': '3d(yz)', 
        '3d_x2-y2': '3d(x^2-y^2)', '3d_z2': '3d(z^2)',
        '4s': '4s', '4px': '4px', '4py': '4py', '4pz': '4pz',
        '4d_xy': '4d(xy)', '4d_xz': '4d(xz)', '4d_yz': '4d(yz)', 
        '4d_x2-y2': '4d(x^2-y^2)', '4d_z2': '4d(z^2)',
        '4f_xyz': '4f(xyz)', '4f_xz2': '4f(xz^2)', '4f_yz2': '4f(yz^2)', 
        '4f_z3': '4f(z^3)', '4f_z(x2-y2)': '4f(z(x^2-y^2))', 
        '4f_x(x2-3y2)': '4f(x(x^2-3y^2))', '4f_y(3x2-y2)': '4f(y(3x^2-y^2))',
        // n=5 轨道
        '5s': '5s', '5px': '5px', '5py': '5py', '5pz': '5pz',
        '5d_xy': '5d(xy)', '5d_xz': '5d(xz)', '5d_yz': '5d(yz)',
        '5d_x2-y2': '5d(x^2-y^2)', '5d_z2': '5d(z^2)',
        '5f_xyz': '5f(xyz)', '5f_xz2': '5f(xz^2)', '5f_yz2': '5f(yz^2)',
        '5f_z3': '5f(z^3)', '5f_z(x2-y2)': '5f(z(x^2-y^2))',
        '5f_x(x2-3y2)': '5f(x(x^2-3y^2))', '5f_y(3x2-y2)': '5f(y(3x^2-y^2))',
        // 5g 轨道
        '5g_z4': '5g(z^4)', '5g_z3x': '5g(z^3·x)', '5g_z3y': '5g(z^3·y)',
        '5g_z2xy': '5g(z^2·xy)', '5g_z2(x2-y2)': '5g(z^2(x^2-y^2))',
        '5g_zx(x2-3y2)': '5g(zx(x^2-3y^2))', '5g_zy(3x2-y2)': '5g(zy(3x^2-y^2))',
        '5g_xy(x2-y2)': '5g(xy(x^2-y^2))', '5g_x4-6x2y2+y4': '5g(x^4-6x^2y^2+y^4)'
    }
};

// UI 元素引用
window.ElectronCloud.ui = {
    orbitalSelect: document.getElementById('orbital-select'),
    axesSizeRange: document.getElementById('axes-size-range'),
    angular3dToggle: document.getElementById('angular-3d-toggle'),
    pointCountSpan: document.getElementById('point-count'),
    startButton: document.getElementById('start-button'),
    maxPointsSelect: document.getElementById('max-points-select'),
    speedRange: document.getElementById('speed-range'),
    centerLock: document.getElementById('center-lock'),
    sizeRange: document.getElementById('size-range'),
    opacityRange: document.getElementById('opacity-range'),
    pauseButton: document.getElementById('pause-button'),
    clearButton: document.getElementById('clear-button'),
    plotTypeSelect: document.getElementById('plot-type-select'),
    multiselectToggle: document.getElementById('multiselect-toggle'),
    compareToggle: document.getElementById('compare-toggle'),
    clearAllSelectionsBtn: document.getElementById('clear-all-selections')
};

// 初始化函数
window.ElectronCloud.init = function() {
    // 初始化 Web Worker 池（用于并行采样加速）
    if (window.ElectronCloud.Sampling.initWorkers) {
        window.ElectronCloud.Sampling.initWorkers();
    }
    
    // 初始化场景
    window.ElectronCloud.Scene.init();
    
    // 初始化UI
    window.ElectronCloud.UI.init();
    
    // 初始化轨道管理
    window.ElectronCloud.Orbital.init();
    
    // 启动动画循环
    window.ElectronCloud.Scene.animate();
    
    // 初始化角向分布3D开关状态
    window.ElectronCloud.UI.updateAngular3DToggleState();
    
    // 确保初始状态下轨道选择正常工作
    const orbitalSelect = window.ElectronCloud.ui.orbitalSelect;
    if (orbitalSelect && orbitalSelect.selectedOptions.length === 0) {
        orbitalSelect.options[0].selected = true;
        const changeEvent = new Event('change', { bubbles: true });
        orbitalSelect.dispatchEvent(changeEvent);
    }
    
    console.log('ElectronCloud 核心模块初始化完成');
};

// 全局状态管理
window.ElectronCloud.setState = function(key, value) {
    window.ElectronCloud.state[key] = value;
};

window.ElectronCloud.getState = function(key) {
    return window.ElectronCloud.state[key];
};

// 重置所有状态（完全重置，用于"重置"按钮）
window.ElectronCloud.resetState = function() {
    const state = window.ElectronCloud.state;
    
    // 重置 Worker 采样会话，使旧的 Worker 结果无效
    if (window.ElectronCloud.Sampling && window.ElectronCloud.Sampling.resetSamplingSession) {
        window.ElectronCloud.Sampling.resetSamplingSession();
    }
    
    // 停止动画
    if (state.animationFrameId) {
        cancelAnimationFrame(state.animationFrameId);
        state.animationFrameId = null;
    }
    
    // 重置采样相关状态
    state.pointCount = 0;
    state.isDrawing = false;
    state.samplingCompleted = false;
    state.farthestDistance = 0;
    state.samplingBoundary = 40;
    state.angularUpdated = false;
    state.axesScaleFactor = 0;
    state.roundRobinIndex = 0; // 重置轮流采样索引
    state.samplingStartTime = null; // 重置采样开始时间
    state.lastChartUpdateTime = 0; // 重置图表更新时间
    
    // 清空数据数组
    state.radialSamples = [];
    state.angularSamples = [];
    
    // 重置轨道状态
    state.orbitalVisibility = {};
    state.renderingCompleted = false;
    state.orbitalPointsMap = {};
    state.originalPositions = null;
    state.orbitalSamplesMap = {};
    state.pointOrbitalIndices = null;
    
    // 重置星空闪烁的基础颜色
    state.baseColors = null;
    state.baseColorsCount = 0;
    
    // 重置闪烁模式缓存（扩散模式和波浪模式）
    state.diffuseDensitiesComputed = false;
    state.waveRanksComputed = false;
    state.waveRanksPointCount = 0;
    state.waveRanks = null;
    
    // 重置自动旋转状态
    if (state.autoRotate) {
        state.autoRotate.enabled = false;
        state.autoRotate.totalAngle = 0;
    }
    
    // 如果正在录制，停止录制
    if (state.isRecordingRotation) {
        state.isRecordingRotation = false;
        if (window.ElectronCloud.UI && window.ElectronCloud.UI.stopRotationRecording) {
            window.ElectronCloud.UI.stopRotationRecording();
        }
    }
    
    // 重置图表数据
    state.backgroundChartData = {
        radial: null,
        angular: null
    };
    
    // 更新UI显示
    if (window.ElectronCloud.ui.pointCountSpan) {
        window.ElectronCloud.ui.pointCountSpan.textContent = '0';
    }
    if (window.ElectronCloud.ui.axesSizeRange) {
        window.ElectronCloud.ui.axesSizeRange.value = '0';
        // 实际重置坐标轴显示
        window.ElectronCloud.UI.onAxesSizeChange({ target: window.ElectronCloud.ui.axesSizeRange });
    }
    
    // 确保坐标系在重置后隐藏
    if (state.customAxes) {
        state.customAxes.visible = false;
    }
    
    console.log('ElectronCloud 状态已完全重置');
};

// 重置采样状态（保留用户设置，用于"启动"按钮）
window.ElectronCloud.resetSamplingState = function() {
    const state = window.ElectronCloud.state;
    
    // 重置 Worker 采样会话，使旧的 Worker 结果无效
    if (window.ElectronCloud.Sampling && window.ElectronCloud.Sampling.resetSamplingSession) {
        window.ElectronCloud.Sampling.resetSamplingSession();
    }
    
    // 停止动画
    if (state.animationFrameId) {
        cancelAnimationFrame(state.animationFrameId);
        state.animationFrameId = null;
    }
    
    // 重置采样相关状态，但保留用户设置（如axesScaleFactor）
    state.pointCount = 0;
    state.isDrawing = false;
    state.samplingCompleted = false;
    state.farthestDistance = 0;
    state.samplingBoundary = 40;
    state.angularUpdated = false;
    state.roundRobinIndex = 0; // 重置轮流采样索引
    state.samplingStartTime = null; // 重置采样开始时间
    state.lastChartUpdateTime = 0; // 重置图表更新时间
    // 从滑动条同步axesScaleFactor，保留用户设置
    if (window.ElectronCloud.ui.axesSizeRange) {
        const sliderValue = parseInt(window.ElectronCloud.ui.axesSizeRange.value, 10);
        state.axesScaleFactor = sliderValue / 100;
    }
    
    // 清空数据数组
    state.radialSamples = [];
    state.angularSamples = [];
    
    // 重置轨道状态
    state.orbitalVisibility = {};
    state.renderingCompleted = false;
    state.orbitalPointsMap = {};
    state.originalPositions = null;
    state.orbitalSamplesMap = {};
    state.pointOrbitalIndices = null;
    
    // 重置星空闪烁的基础颜色
    state.baseColors = null;
    state.baseColorsCount = 0;
    
    // 重置闪烁模式缓存（扩散模式和波浪模式）
    state.diffuseDensitiesComputed = false;
    state.waveRanksComputed = false;
    state.waveRanksPointCount = 0;
    state.waveRanks = null;
    
    // 重置图表数据
    state.backgroundChartData = {
        radial: null,
        angular: null
    };
    
    // 更新UI显示
    if (window.ElectronCloud.ui.pointCountSpan) {
        window.ElectronCloud.ui.pointCountSpan.textContent = '0';
    }
    // 不重置坐标轴滑动条值，保留用户设置
    
    // 由于farthestDistance重置为0，需要隐藏坐标系直到有新的采样点
    if (state.customAxes) {
        state.customAxes.visible = false;
    }
    
    console.log('ElectronCloud 采样状态已重置，保留用户设置');
};
