// 轨道管理和数据处理模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Orbital = {};

// 初始化轨道管理
window.ElectronCloud.Orbital.init = function() {
    console.log('轨道管理模块初始化完成');
};

// 开始绘制轨道
window.ElectronCloud.Orbital.startDrawing = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    // 读取轨道参数
    const selected = Array.from(ui.orbitalSelect.selectedOptions || []).map(o => o.value);
    
    // 【关键检查】在多选或比照模式下，必须选择至少一个轨道才能开始渲染
    const isMultiselectMode = ui.multiselectToggle && ui.multiselectToggle.checked;
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    
    if ((isMultiselectMode || isCompareMode) && selected.length === 0) {
        alert('请先选择至少一个轨道再启动渲染');
        return;
    }
    
    // 停止当前绘制
    state.isDrawing = false;
    
    if (state.animationFrameId) {
        cancelAnimationFrame(state.animationFrameId);
    }

    // 清除旧的点云
    window.ElectronCloud.Scene.clearPoints();
    
    // 重置采样状态（保留用户设置）
    window.ElectronCloud.resetSamplingState();
    
    // 重置角向更新标记
    state.angularUpdated = false;
    
    // 设置当前轨道
    state.currentOrbitals = selected.length ? selected : [ui.orbitalSelect.value || '1s'];
    state.currentOrbital = state.currentOrbitals[0];
    
    // 动态调整采样边界：根据所选轨道的最大 n 值
    let maxN = 1;
    for (const key of state.currentOrbitals) {
        const p = Hydrogen.orbitalParamsFromKey(key);
        if (p && p.n > maxN) maxN = p.n;
    }
    // 初始边界使用高效公式：4n² + 2n，覆盖约95-99%概率区域
    // 剩余尾部由动态扩展机制自动捕获（采样到远处点时边界自动扩大）
    // n=1 -> 6, n=2 -> 20, n=3 -> 42, n=4 -> 72, n=5 -> 110
    // 对于高 n 值轨道，使用重要性采样可以更高效地覆盖弥散区域
    state.samplingBoundary = Math.max(12, 4 * maxN * maxN + 2 * maxN);
    console.log(`根据最大主量子数 n=${maxN} 调整采样边界为: ${state.samplingBoundary}`);
    
    const maxPointsValue = ui.maxPointsSelect.value;
    state.MAX_POINTS = parseInt(maxPointsValue, 10);

    // 创建新的点集合
    const bufferSize = state.MAX_POINTS;
    const geometry = new THREE.BufferGeometry();
    const positions = new Float32Array(bufferSize * 3);
    const colors = new Float32Array(bufferSize * 3); // 预创建颜色数组
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3)); // 添加颜色属性

    // 创建材质
    const sprite = window.ElectronCloud.Scene.generateCircleSprite();
    
    // 从亮度滑条计算初始透明度 (0-80 -> 0.05-1.0)
    const brightnessValue = parseInt(ui.opacityRange.value, 10);
    const initialOpacity = brightnessValue <= 80 
        ? 0.05 + (brightnessValue / 80) * 0.95 
        : 1.0;
    
    const material = new THREE.PointsMaterial({
        map: sprite,
        size: window.ElectronCloud.UI.getPointSize(),
        transparent: true,
        depthWrite: false,
        opacity: initialOpacity,
        vertexColors: true
    });

    state.points = new THREE.Points(geometry, material);
    state.scene.add(state.points);
    
    // 【根本修复】统一重置所有场景对象的旋转状态
    // 这是解决"坐标歪斜"问题的关键步骤
    // 确保新的点云与坐标轴、角向叠加层三者旋转同步（都归零）
    if (window.ElectronCloud.Scene.resetAllSceneObjectsRotation) {
        window.ElectronCloud.Scene.resetAllSceneObjectsRotation();
    }
    
    // 【修改】渲染开始时，只重置totalAngle，保留用户预设的enabled状态
    // 这样用户可以在渲染前点击"自动旋转"预设，渲染后自动应用
    if (state.autoRotate) {
        state.autoRotate.totalAngle = 0;
        // 不重置enabled，保留用户预设的状态
    }
    // 录制按钮在自动旋转未启用时仍然禁用
    const recordBtn = document.getElementById('record-rotation-btn');
    if (recordBtn) recordBtn.disabled = true;
    
    // 重置锁定视角的UI状态
    const lockButtons = document.querySelectorAll('.lock-axis-btn');
    lockButtons.forEach(btn => {
        btn.classList.remove('active', 'green');
    });
    
    // 注意：resetAllSceneObjectsRotation 已经重置了 customAxes 和 angularOverlay
    // 新创建的 points 不需要重置，因为刚创建时旋转就是 (0,0,0)

    // 开始绘制
    state.isDrawing = true;
    state.samplingStartTime = performance.now(); // 记录采样开始时间
    state.lastChartUpdateTime = 0; // 重置图表更新时间，确保新轨道能正常刷新
    
    // 更新坐标轴滑动条状态（渲染中置灰）
    window.ElectronCloud.UI.updateAxesSizeRangeState();
    
    console.log('开始采样，初始化实时图表');
    
    // 预计算理论曲线数据（不需要采样数据）
    window.ElectronCloud.Orbital.initTheoryData();
    
    // 确保数据面板可见（如果用户之前滑出收起了）
    const dataPanel = document.getElementById('data-panel');
    if (dataPanel) {
        dataPanel.classList.remove('collapsed');
    }
    
    // 更新角向分布3D开关状态
    window.ElectronCloud.UI.updateAngular3DToggleState();
    
    // 更新清空按钮状态
    window.ElectronCloud.UI.updateClearAllSelectionsState();

    // 更新全屏按钮状态
    window.ElectronCloud.UI.updateFullscreenBtnState();
    
    // 禁用手势按钮
    window.ElectronCloud.UI.disableGestureButton();
    
    // 更新角向分布叠加
    const angular3dToggle = ui.angular3dToggle;
    if (angular3dToggle && angular3dToggle.checked) {
        window.ElectronCloud.Visualization.updateAngularOverlay();
    }
    
    // 如果用户已设置坐标轴比例系数，立即应用
    if (state.axesScaleFactor > 0 && ui.axesSizeRange) {
        window.ElectronCloud.UI.onAxesSizeChange({ target: ui.axesSizeRange });
    }
    
    window.ElectronCloud.Scene.animate();
};

// 清除绘制
window.ElectronCloud.Orbital.clearDrawing = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    // 停止绘制
    state.isDrawing = false;
    
    // 更新坐标轴滑动条状态（可修改）
    window.ElectronCloud.UI.updateAxesSizeRangeState();

    if (state.animationFrameId) {
        cancelAnimationFrame(state.animationFrameId);
    }

    // 【根本修复】统一重置所有场景对象的旋转状态
    // 确保重置后坐标轴、点云（即将清除）、角向叠加层旋转归零
    if (window.ElectronCloud.Scene.resetAllSceneObjectsRotation) {
        window.ElectronCloud.Scene.resetAllSceneObjectsRotation();
    }

    // 清除点云
    window.ElectronCloud.Scene.clearPoints();
    
    // 重置状态
    window.ElectronCloud.resetState();
    
    // 重置所有闪烁模式的密度缓存
    state.diffuseDensitiesComputed = false;
    state.waveRanksComputed = false;
    state.waveRanksPointCount = 0;
    state.waveRanks = null;
    
    // 重置轨道选项UI状态 - 清除所有可见性相关的类和样式
    const options = ui.orbitalSelect.querySelectorAll('option');
    options.forEach(option => {
        option.classList.remove('orbital-hidden', 'orbital-visible');
        option.style.opacity = '';
        option.style.textDecoration = '';
        option.title = '';
    });
    
    // 【修复】同时清除自定义轨道列表的隐藏样式
    const customList = document.getElementById('custom-orbital-list');
    if (customList) {
        Array.from(customList.children).forEach(item => {
            item.classList.remove('orbital-hidden', 'compare-a', 'compare-b', 'compare-c', 'active');
        });
    }
    
    // 刷新选择框样式，清除比照模式的颜色状态
    window.ElectronCloud.UI.refreshSelectStyles();
    
    // 如果当前在比照模式或多选模式，更新计数显示
    if ((ui.multiselectToggle && ui.multiselectToggle.checked) || (ui.compareToggle && ui.compareToggle.checked)) {
        window.ElectronCloud.UI.updateSelectionCount();
    }
    
    // 移除渲染完成状态的CSS类
    const controlPanel = document.getElementById('control-panel');
    if (controlPanel) {
        controlPanel.classList.remove('rendering-completed');
    }
    
    // 更新清空按钮状态（重置后重新启用）
    window.ElectronCloud.UI.updateClearAllSelectionsState();

    // 更新全屏按钮状态
    window.ElectronCloud.UI.updateFullscreenBtnState();
    
    console.log('清除绘制');
    
    // 重置时不改变面板的展开/收起状态，用户保留当前状态

    // 清除角向分布叠加
    window.ElectronCloud.Scene.clearAngularOverlay();

    // 重启动画循环以保持控制器工作
    window.ElectronCloud.Scene.animate();
    
    // 更新角向分布3D开关状态
    window.ElectronCloud.UI.updateAngular3DToggleState();
    
    // 重置自动旋转 UI 状态
    const autoRotateToggle = document.getElementById('auto-rotate-toggle');
    const rotationFeatureBox = document.getElementById('rotation-feature-box');
    if (autoRotateToggle) {
        autoRotateToggle.checked = false;
        autoRotateToggle.dispatchEvent(new Event('change', { bubbles: true }));
    }
    if (rotationFeatureBox) {
        rotationFeatureBox.classList.remove('active');
    }
    
    // 【新增】重置闪烁模式 UI 状态
    const heartbeatToggle = document.getElementById('heartbeat-toggle');
    const flickerFeatureBox = document.getElementById('flicker-feature-box');
    if (heartbeatToggle && heartbeatToggle.checked) {
        heartbeatToggle.checked = false;
        heartbeatToggle.dispatchEvent(new Event('change', { bubbles: true }));
    }
    if (flickerFeatureBox) {
        flickerFeatureBox.classList.remove('active');
    }
    
    // 重置锁定视角的UI状态
    const lockButtons = document.querySelectorAll('.lock-axis-btn');
    lockButtons.forEach(btn => {
        btn.classList.remove('active', 'green');
    });
    // 更新自动旋转按钮状态（清除后保持可用，允许预设）
    window.ElectronCloud.UI.updateAutoRotateButtonState();
};

// 更新轨道可见性
window.ElectronCloud.Orbital.updateOrbitalVisibility = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    // 只有渲染完成后才能切换可见性
    if (!state.points || !state.renderingCompleted || !ui.compareToggle.checked) return;
    
    const positions = state.points.geometry.attributes.position.array;
    
    // 如果还没有备份原始位置，先备份
    if (!state.originalPositions) {
        state.originalPositions = new Float32Array(positions);
    }
    
    // 【性能优化】预先构建点索引到轨道的映射，避免 O(n²) 复杂度
    // 使用 pointOrbitalIndices 数组（在采样时已经记录）
    const orbitalKeys = state.currentOrbitals || [];
    
    // 遍历所有点，决定是否显示
    for (let i = 0; i < state.pointCount; i++) {
        const orbitalIndex = state.pointOrbitalIndices ? state.pointOrbitalIndices[i] : 0;
        const orbitalKey = orbitalKeys[orbitalIndex];
        const isVisible = orbitalKey ? (state.orbitalVisibility[orbitalKey] !== false) : true;
        
        const posIdx = i * 3;
        if (isVisible) {
            // 显示：恢复原始位置
            positions[posIdx] = state.originalPositions[posIdx];
            positions[posIdx + 1] = state.originalPositions[posIdx + 1];
            positions[posIdx + 2] = state.originalPositions[posIdx + 2];
        } else {
            // 隐藏：移动到视野外
            positions[posIdx] = 10000;
            positions[posIdx + 1] = 10000;
            positions[posIdx + 2] = 10000;
        }
    }
    
    state.points.geometry.attributes.position.needsUpdate = true;
};

// 切换轨道可见性
window.ElectronCloud.Orbital.toggleOrbitalVisibility = function(orbitalKey) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    console.log(`toggleOrbitalVisibility called with: ${orbitalKey}`);
    console.log(`renderingCompleted: ${state.renderingCompleted}, compareToggle.checked: ${ui.compareToggle?.checked}`);
    console.log(`currentOrbitals: ${JSON.stringify(state.currentOrbitals)}`);
    
    // 只有渲染完成后才能切换可见性
    if (!state.renderingCompleted || !ui.compareToggle.checked) {
        console.log('渲染未完成或不在比照模式，返回');
        return;
    }
    
    // 只有已渲染的轨道才能被切换可见性
    if (!state.currentOrbitals.includes(orbitalKey)) {
        console.log(`轨道 ${orbitalKey} 未被渲染，无法切换可见性。currentOrbitals: ${state.currentOrbitals.join(', ')}`);
        return;
    }
    
    console.log(`开始切换轨道 ${orbitalKey} 的可见性`);
    
    // 确保 orbitalVisibility 中有这个轨道的记录
    if (state.orbitalVisibility[orbitalKey] === undefined) {
        state.orbitalVisibility[orbitalKey] = true;
    }
    
    state.orbitalVisibility[orbitalKey] = !state.orbitalVisibility[orbitalKey];
    console.log(`轨道 ${orbitalKey} 可见性设为: ${state.orbitalVisibility[orbitalKey]}`);
    
    window.ElectronCloud.Orbital.updateOrbitalVisibility();
    
    // 更新UI显示 - 使用删除线和透明度来表示隐藏状态，但保持选中状态不变
    const option = ui.orbitalSelect.querySelector(`option[value="${orbitalKey}"]`);
    if (option) {
        if (state.orbitalVisibility[orbitalKey]) {
            option.classList.remove('orbital-hidden');
            option.classList.add('orbital-visible');
            option.title = '点击隐藏此轨道';
        } else {
            option.classList.add('orbital-hidden');
            option.classList.remove('orbital-visible');
            option.title = '点击显示此轨道';
        }
    }
    
    // 同步更新图表中对应曲线的可见性
    window.ElectronCloud.Orbital.syncChartVisibility(orbitalKey);
};

// 同步图表曲线的可见性
window.ElectronCloud.Orbital.syncChartVisibility = function(orbitalKey) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    
    // 获取图表实例
    const chart = window.DataPanel?.state?.chart;
    if (!chart || !chart.data || !chart.data.datasets) {
        console.log('图表不存在或没有数据集');
        return;
    }
    
    // 使用统一的轨道显示名称映射
    const displayName = constants.orbitalDisplayNames[orbitalKey] || orbitalKey;
    const isVisible = state.orbitalVisibility[orbitalKey] !== false;
    
    console.log(`同步图表曲线: orbitalKey=${orbitalKey}, displayName=${displayName}, isVisible=${isVisible}`);
    console.log(`图表数据集数量: ${chart.data.datasets.length}`);
    console.log(`数据集标签: ${chart.data.datasets.map(ds => ds.label).join(', ')}`);
    
    // 查找并更新对应曲线的hidden状态
    let found = false;
    for (let i = 0; i < chart.data.datasets.length; i++) {
        const dataset = chart.data.datasets[i];
        if (dataset.label === displayName) {
            dataset.hidden = !isVisible;
            found = true;
            console.log(`找到匹配的数据集，设置 hidden=${!isVisible}`);
            break;
        }
    }
    
    if (!found) {
        console.log(`未找到匹配的数据集: ${displayName}`);
    }
    
    // 更新图表
    chart.update('none');
};

// 初始化理论曲线数据（在采样开始前调用，显示参考线）
window.ElectronCloud.Orbital.initTheoryData = function() {
    const state = window.ElectronCloud.state;
    
    if (!window.Hydrogen) return;
    
    try {
        // 【关键修复】直接使用 state.currentOrbitals，不依赖 UI 状态判断
        // state.currentOrbitals 在 startDrawing 时已经被正确设置为所有选中的轨道
        const orbitals = state.currentOrbitals && state.currentOrbitals.length > 0 
            ? state.currentOrbitals 
            : [state.currentOrbital || '1s'];
        
        // 获取所有轨道的参数
        const paramsList = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
        if (paramsList.length === 0) return;
        
        console.log('initTheoryData: 计算理论曲线，轨道列表:', orbitals);
        
        const angularBins = 90;
        
        // 使用最大的 n 来计算预估最大半径
        const maxN = Math.max(...paramsList.map(p => p.n));
        const estimatedRmax = Hydrogen.recommendRmax(maxN);
        const radialBins = 240;
        
        // 创建空的径向直方图（零值）
        const radialEdges = new Float32Array(radialBins + 1);
        const radialCounts = new Float32Array(radialBins);
        const dr = estimatedRmax / radialBins;
        for (let i = 0; i <= radialBins; i++) radialEdges[i] = i * dr;
        
        const radialCenters = new Array(radialBins);
        const radialTheoryValues = new Array(radialBins);
        const radialWave = new Array(radialBins);
        
        // 【关键修复】多选模式下，理论曲线是各轨道概率密度的加和
        for (let i = 0; i < radialBins; i++) {
            radialCenters[i] = 0.5 * (radialEdges[i] + radialEdges[i + 1]);
            
            // 对所有轨道的概率密度求和
            let sumPDF = 0;
            let sumWave = 0;
            for (const params of paramsList) {
                sumPDF += Hydrogen.radialPDF(params.n, params.l, radialCenters[i]);
                sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, radialCenters[i]));
            }
            // 多轨道时归一化（除以轨道数，因为采样时各轨道等概率）
            radialTheoryValues[i] = sumPDF / paramsList.length;
            radialWave[i] = sumWave / paramsList.length;
        }
        
        state.backgroundChartData.radial = {
            hist: { edges: radialEdges, counts: radialCounts, dr: dr, rmax: estimatedRmax },
            theory: { centers: radialCenters, values: radialTheoryValues, wave: radialWave }
        };
        
        // 创建空的角向直方图（零值）
        const angularEdges = new Float32Array(angularBins + 1);
        const angularCounts = new Float32Array(angularBins);
        const dth = Math.PI / angularBins;
        for (let i = 0; i <= angularBins; i++) angularEdges[i] = i * dth;
        
        const angularCenters = new Array(angularBins);
        const angularTheoryValues = new Array(angularBins);
        
        // 【关键修复】多选模式下，角向理论曲线也是各轨道的加和
        for (let i = 0; i < angularBins; i++) {
            angularCenters[i] = 0.5 * (angularEdges[i] + angularEdges[i + 1]);
            
            let sumAngularPDF = 0;
            for (const params of paramsList) {
                sumAngularPDF += Hydrogen.angularPDF_Theta(params.angKey.l, params.angKey.m, angularCenters[i]);
            }
            angularTheoryValues[i] = sumAngularPDF / paramsList.length;
        }
        
        state.backgroundChartData.angular = {
            hist: { edges: angularEdges, counts: angularCounts, 'dθ': dth },
            theory: { centers: angularCenters, values: angularTheoryValues }
        };
        
        console.log('理论曲线数据初始化完成，轨道数:', paramsList.length);
        
    } catch (error) {
        console.error('初始化理论数据失败:', error);
    }
};

// 后台准备图表数据
window.ElectronCloud.Orbital.updateBackgroundChartData = function() {
    const state = window.ElectronCloud.state;
    
    if (!window.Hydrogen || state.radialSamples.length === 0) return;
    
    try {
        const angularBins = 90;
        
        // 【关键修复】直接使用 state.currentOrbitals，不依赖 UI 状态判断
        const orbitals = state.currentOrbitals && state.currentOrbitals.length > 0 
            ? state.currentOrbitals 
            : [state.currentOrbital || '1s'];
        
        // 获取所有轨道的参数
        const paramsList = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
        if (paramsList.length === 0) return;

        // 准备径向数据
        if (state.radialSamples.length > 0) {
            // 【性能修复】使用循环替代Math.max(...array)，避免大数组栈溢出
            let maxSampleDistance = 0;
            for (let i = 0; i < state.radialSamples.length; i++) {
                if (state.radialSamples[i] > maxSampleDistance) {
                    maxSampleDistance = state.radialSamples[i];
                }
            }
            const dynamicRmax = Math.max(1, maxSampleDistance * 1.08);
            
            // 动态调整箱数
            const baseBins = 240;
            const sampleDensity = state.radialSamples.length / Math.max(1, dynamicRmax);
            const adaptiveBins = Math.min(400, Math.max(baseBins, Math.floor(sampleDensity * 0.5)));
            
            const radialHist = Hydrogen.histogramRadialFromSamples(state.radialSamples, adaptiveBins, dynamicRmax, true);
            
            // 横轴中心
            const centers = new Array(adaptiveBins);
            for (let i = 0; i < adaptiveBins; i++) centers[i] = 0.5 * (radialHist.edges[i] + radialHist.edges[i + 1]);
            
            // 计算理论值 - 多轨道时加和后取平均
            const values = new Array(adaptiveBins);
            const wave = new Array(adaptiveBins);
            for (let i = 0; i < adaptiveBins; i++) {
                let sumPDF = 0;
                let sumWave = 0;
                for (const params of paramsList) {
                    sumPDF += Hydrogen.radialPDF(params.n, params.l, centers[i]);
                    sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, centers[i]));
                }
                values[i] = sumPDF / paramsList.length;
                wave[i] = sumWave / paramsList.length;
            }

            state.backgroundChartData.radial = {
                hist: radialHist,
                theory: { centers, values, wave }
            };
        }

        // 准备角向数据
        if (state.angularSamples.length > 0) {
            const angularHist = Hydrogen.histogramThetaFromSamples(state.angularSamples, angularBins, true);
            
            // 【关键修复】计算角向理论值 - 多轨道时加和
            const centers = new Array(angularBins);
            const values = new Array(angularBins);
            for (let i = 0; i < angularBins; i++) {
                centers[i] = 0.5 * (angularHist.edges[i] + angularHist.edges[i+1]);
                
                let sumAngularPDF = 0;
                for (const params of paramsList) {
                    sumAngularPDF += Hydrogen.angularPDF_Theta(params.angKey.l, params.angKey.m, centers[i]);
                }
                values[i] = sumAngularPDF / paramsList.length;
            }
            
            state.backgroundChartData.angular = {
                hist: angularHist,
                theory: { centers, values }
            };
        }
        
        console.log('后台图表数据已更新，径向数据点:', state.radialSamples.length, '轨道数:', paramsList.length);
        
    } catch (error) {
        console.error('后台更新图表数据失败:', error);
    }
};

// 绘制概率图表
window.ElectronCloud.Orbital.drawProbabilityChart = function(final = true) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    if (!window.DataPanel || !window.Hydrogen) return;
    
    console.log('开始绘制概率图表，final:', final);
    
    // 确保数据面板可见（不改变展开/收起状态）
    const panel = document.getElementById('data-panel');
    if (panel) {
        panel.classList.remove('collapsed');
    }

    const type = ui.plotTypeSelect ? ui.plotTypeSelect.value : 'radial';
    
    // 检查是否是对比模式
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    
    if (isCompareMode && Object.keys(state.orbitalSamplesMap).length > 0) {
        console.log('使用对比模式散点图渲染，轨道数量:', Object.keys(state.orbitalSamplesMap).length);
        DataPanel.renderChartCompare(state.orbitalSamplesMap, type);
        return;
    }
    
    // 优先使用后台准备的数据
    if (type === 'radial' && state.backgroundChartData.radial) {
        console.log('使用后台准备的径向图表数据进行渲染');
        DataPanel.renderChartRadial(state.backgroundChartData.radial.hist, state.backgroundChartData.radial.theory);
        return;
    } else if (type === 'angular' && state.backgroundChartData.angular) {
        console.log('使用后台准备的角向图表数据进行渲染');
        DataPanel.renderChartAngular(state.backgroundChartData.angular.hist, state.backgroundChartData.angular.theory);
        return;
    }
    
    // 如果没有后台数据，则重新计算
    console.log('后台数据不可用，重新计算图表数据');
    const angularBins = 90;
    
    // 【关键修复】使用所有选中的轨道计算理论曲线，而不是单个轨道
    const orbitals = state.currentOrbitals && state.currentOrbitals.length > 0 
        ? state.currentOrbitals 
        : [state.currentOrbital || '1s'];
    const paramsList = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
    
    if (paramsList.length === 0) return;

    if (type === 'radial') {
        // 动态调整横轴宽度：仅略宽于最远采样点距离
        // 【性能修复】使用循环替代Math.max(...array)，避免大数组栈溢出
        let maxSampleDistance = 0;
        for (let i = 0; i < state.radialSamples.length; i++) {
            if (state.radialSamples[i] > maxSampleDistance) {
                maxSampleDistance = state.radialSamples[i];
            }
        }
        const dynamicRmax = Math.max(1, maxSampleDistance * 1.08);
        
        // 动态调整箱数：确保足够细致，根据数据范围和采样点数自适应
        const baseBins = 240;
        const sampleDensity = state.radialSamples.length / Math.max(1, dynamicRmax);
        const adaptiveBins = Math.min(400, Math.max(baseBins, Math.floor(sampleDensity * 0.5)));
        
        const hist = Hydrogen.histogramRadialFromSamples(state.radialSamples, adaptiveBins, dynamicRmax, true);
        
        // 横轴中心
        const centers = new Array(adaptiveBins);
        for (let i = 0; i < adaptiveBins; i++) centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);
        
        // 【关键修复】计算理论值 - 多轨道时加和
        const values = new Array(adaptiveBins);
        const wave = new Array(adaptiveBins);
        for (let i = 0; i < adaptiveBins; i++) {
            let sumPDF = 0;
            let sumWave = 0;
            for (const params of paramsList) {
                sumPDF += Hydrogen.radialPDF(params.n, params.l, centers[i]);
                sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, centers[i]));
            }
            values[i] = sumPDF / paramsList.length;
            wave[i] = sumWave / paramsList.length;
        }

        // 绘制直方图和理论曲线
        DataPanel.renderChartRadial(hist, { centers, values, wave });
    } else {
        const hist = Hydrogen.histogramThetaFromSamples(state.angularSamples, angularBins, true);
        
        // 【关键修复】计算角向理论值 - 多轨道时加和
        const centers = new Array(angularBins);
        const values = new Array(angularBins);
        for (let i = 0; i < angularBins; i++) {
            centers[i] = 0.5 * (hist.edges[i] + hist.edges[i+1]);
            
            let sumAngularPDF = 0;
            for (const params of paramsList) {
                sumAngularPDF += Hydrogen.angularPDF_Theta(params.angKey.l, params.angKey.m, centers[i]);
            }
            values[i] = sumAngularPDF / paramsList.length;
        }
        
        DataPanel.renderChartAngular(hist, { centers, values });
    }
};
