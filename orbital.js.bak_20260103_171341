// 轨道管理和数据处理模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Orbital = {};

// 初始化轨道管理
window.ElectronCloud.Orbital.init = function () {
    console.log('轨道管理模块初始化完成');
};

// 开始绘制轨道
window.ElectronCloud.Orbital.startDrawing = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    // 读取轨道参数
    const selected = Array.from(ui.orbitalSelect.selectedOptions || []).map(o => o.value);

    // 【修改】所有模式下都必须选择至少一个轨道才能开始渲染
    // 启动时不再有默认轨道选择，与多选/比照模式保持一致
    if (selected.length === 0) {
        alert('请先选择至少一个轨道再启动渲染');
        return;
    }

    // 杂化模式下验证轨道组合是否为支持的标准构型
    if (state.isHybridMode && !state.hybridizationValid) {
        alert('当前轨道组合不支持杂化计算\n\n仅支持：sp、sp²、sp³、sp³d、sp³d²');
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

    // 设置当前轨道（不再有 fallback，因为上面已经检查过 selected.length）
    state.currentOrbitals = selected;
    state.currentOrbital = selected[0];

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

    // 【新增】根据预估轨道半径自动缩放相机
    // 目标：使轨道轮廓最长处占画面高度的约 2/3
    if (Hydrogen.estimateOrbitalRadius95 && state.camera) {
        let maxRadius = 15; // 默认 1s 轨道
        const atomType = state.currentAtom || 'H';

        for (const key of state.currentOrbitals) {
            const r95 = Hydrogen.estimateOrbitalRadius95(atomType, key);
            if (r95 > maxRadius) maxRadius = r95;
        }

        // 相机距离计算：使轨道占画面高度 2/3
        // FOV = 75°，半FOV = 37.5°
        // 占画面 2/3 意味着轨道半径对应的视角是 (2/3) × 37.5° ≈ 25°
        // distance = radius / tan(25°) ≈ radius × 2.14
        // 但用户反馈太远，改为 radius × 1.2 使轨道更大
        const targetDistance = maxRadius * 0.3;

        // 限制在合理范围内
        const newDistance = Math.max(10, Math.min(500, targetDistance));

        // 保持相机方向不变，只改变距离
        const currentPos = state.camera.position;
        const currentDist = currentPos.length();
        if (currentDist > 0) {
            const scale = newDistance / currentDist;
            state.camera.position.multiplyScalar(scale);
        } else {
            state.camera.position.set(0, 0, newDistance);
        }

        // 更新控制器（如果有）
        if (state.controls) {
            state.controls.update();
        }

        console.log(`自动缩放相机：预估轨道半径=${maxRadius.toFixed(1)}, 相机距离=${newDistance.toFixed(1)}`);
    }

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

    // 【先验缩放】根据轨道尺寸设置辉光半径（渲染前应用）
    // 以氢原子1s轨道为基准，多轨道时使用最小的轨道半径
    if (state.bloomPass && window.ElectronCloud.UI.getBloomRadius) {
        const scaledBloomRadius = window.ElectronCloud.UI.getBloomRadius();
        state.bloomPass.radius = scaledBloomRadius;
        console.log(`先验辉光缩放：Bloom半径设置为 ${scaledBloomRadius.toFixed(3)}`);
    }

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

    // 更新模式切换栏状态（渲染中禁用）
    window.ElectronCloud.UI.updateModeSwitcherState();

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

    // 【新增】采样过程中禁用轨道轮廓开关
    const hybridContourToggle = document.getElementById('hybrid-contour-toggle');
    if (hybridContourToggle) {
        hybridContourToggle.disabled = true;
    }

    // 【新增】采样过程中禁用相位显示开关
    const phaseToggle = document.getElementById('phase-toggle');
    if (phaseToggle) {
        phaseToggle.disabled = true;
        const phaseBox = document.getElementById('phase-box');
        if (phaseBox) {
            phaseBox.classList.add('disabled');
            phaseBox.title = '采样过程中不可切换';
        }
    }

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
window.ElectronCloud.Orbital.clearDrawing = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    // 停止绘制
    state.isDrawing = false;

    // 更新坐标轴滑动条状态（可修改）
    window.ElectronCloud.UI.updateAxesSizeRangeState();

    // 更新模式切换栏状态（可切换）
    window.ElectronCloud.UI.updateModeSwitcherState();

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

    // 【新增】清除绘制时，重新启用轨道轮廓开关
    const hybridContourToggle = document.getElementById('hybrid-contour-toggle');
    if (hybridContourToggle) {
        hybridContourToggle.disabled = false;
    }

    // 更新清空按钮状态（重置后重新启用）
    window.ElectronCloud.UI.updateClearAllSelectionsState();

    // 更新全屏按钮状态
    window.ElectronCloud.UI.updateFullscreenBtnState();

    // 【修复】清除绘制时，真正销毁普通模式的等值面缓存
    const contour3dToggle = document.getElementById('contour-3d-toggle');
    if (contour3dToggle) {
        contour3dToggle.checked = false;
    }
    if (state.contourOverlay) {
        state.scene.remove(state.contourOverlay);
        state.contourOverlay.traverse((child) => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
        state.contourOverlay = null;
    }

    // 【修复】清除绘制时，真正销毁杂化轨道轮廓缓存（因为轨道已改变）
    // 不能只隐藏，必须销毁 geometry 和 material 释放内存
    if (state.hybridContourOverlays && state.hybridContourOverlays.length > 0) {
        for (const overlay of state.hybridContourOverlays) {
            if (!overlay) continue;
            state.scene.remove(overlay);
            overlay.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
        }
        state.hybridContourOverlays = null; // 清空缓存
    }

    // 重置开关状态
    if (hybridContourToggle) {
        hybridContourToggle.checked = false;
        hybridContourToggle.disabled = false;
    }

    // 【新增】清除绘制时，重新启用相位显示开关
    const phaseToggle = document.getElementById('phase-toggle');
    if (phaseToggle) {
        phaseToggle.disabled = false;
        const phaseBox = document.getElementById('phase-box');
        if (phaseBox) {
            phaseBox.classList.remove('disabled');
            phaseBox.title = '';
        }
    }

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
window.ElectronCloud.Orbital.updateOrbitalVisibility = function () {
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
window.ElectronCloud.Orbital.toggleOrbitalVisibility = function (orbitalKey) {
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
window.ElectronCloud.Orbital.syncChartVisibility = function (orbitalKey) {
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
window.ElectronCloud.Orbital.initTheoryData = function () {
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

        const angularBins = 180;

        // 使用最大的 n 来计算预估最大半径
        const maxN = Math.max(...paramsList.map(p => p.n));
        const estimatedRmax = Hydrogen.recommendRmax(maxN);
        const radialBins = 480;

        // 创建空的径向直方图（零值）
        const radialEdges = new Float32Array(radialBins + 1);
        const radialCounts = new Float32Array(radialBins);
        const dr = estimatedRmax / radialBins;
        for (let i = 0; i <= radialBins; i++) radialEdges[i] = i * dr;

        const radialCenters = new Array(radialBins);
        const radialTheoryValues = new Array(radialBins);
        const radialWave = new Array(radialBins);

        // 【关键修复】检测杂化模式，使用正确的物理公式
        const isHybridMode = state.isHybridMode === true;
        const atomType = state.currentAtom || 'H'; // 支持非氢原子
        const atomZ = (window.SlaterBasis && window.SlaterBasis[atomType]) ? window.SlaterBasis[atomType].Z : 1;

        for (let i = 0; i < radialBins; i++) {
            radialCenters[i] = 0.5 * (radialEdges[i] + radialEdges[i + 1]);

            if (isHybridMode) {
                // 【杂化模式】使用专门的杂化径向PDF：P(r) = r² × Σ|cᵢ|²|Rᵢ(r)|²
                radialTheoryValues[i] = Hydrogen.hybridRadialPDF(paramsList, radialCenters[i], atomZ, 1, atomType);
                // 波函数：取加权和
                let sumWave = 0;
                const defaultCoeff = 1.0 / Math.sqrt(paramsList.length);
                for (const params of paramsList) {
                    const coeff = params.coefficient !== undefined ? params.coefficient : defaultCoeff;
                    sumWave += coeff * Math.abs(Hydrogen.radialR(params.n, params.l, radialCenters[i], atomZ, 1, atomType));
                }
                radialWave[i] = sumWave;
            } else {
                // 【普通模式】各轨道PDF的平均
                let sumPDF = 0;
                let sumWave = 0;
                for (const params of paramsList) {
                    sumPDF += Hydrogen.radialPDF(params.n, params.l, radialCenters[i], atomZ, 1, atomType);
                    sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, radialCenters[i], atomZ, 1, atomType));
                }
                radialTheoryValues[i] = sumPDF / paramsList.length;
                radialWave[i] = sumWave / paramsList.length;
            }
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

        // 创建φ角向直方图（理论曲线）
        const phiBins = 180;
        const phiEdges = new Float32Array(phiBins + 1);
        const phiCounts = new Float32Array(phiBins);
        const dPhi = (2 * Math.PI) / phiBins;
        for (let i = 0; i <= phiBins; i++) phiEdges[i] = i * dPhi;

        const phiCenters = new Array(phiBins);
        const phiTheoryValues = new Array(phiBins);

        for (let i = 0; i < phiBins; i++) {
            phiCenters[i] = 0.5 * (phiEdges[i] + phiEdges[i + 1]);
            let sumPhiPDF = 0;
            for (const params of paramsList) {
                sumPhiPDF += Hydrogen.angularPDF_Phi(params.angKey.m, params.angKey.t, phiCenters[i]);
            }
            phiTheoryValues[i] = sumPhiPDF / paramsList.length;
        }

        state.backgroundChartData.phi = {
            hist: { edges: phiEdges, counts: phiCounts, 'dφ': dPhi },
            theory: { centers: phiCenters, values: phiTheoryValues }
        };

        // ==================== 能量数据预计算（简化版：仅使用 ε·P(r)）====================
        try {
            const numBins = radialBins;
            const theoryE = new Float32Array(numBins);
            const theoryDEdr = new Float32Array(numBins);
            let avgEpsilon = -0.5;

            if (window.Hydrogen.calculateCumulativeOrbitalEnergy) {
                for (const p of paramsList) {
                    const pZ = window.SlaterBasis && window.SlaterBasis[atomType]
                        ? window.SlaterBasis[atomType].Z
                        : 1;
                    const res = window.Hydrogen.calculateCumulativeOrbitalEnergy(
                        p.n, p.l, pZ, atomType, radialCenters
                    );
                    if (!res) continue;
                    for (let j = 0; j < numBins; j++) {
                        theoryE[j] += (res.E && res.E.length > j) ? res.E[j] : 0;
                        theoryDEdr[j] += (res.dEdr && res.dEdr.length > j) ? res.dEdr[j] : 0;
                    }
                    if (res.epsilon !== undefined) {
                        avgEpsilon = res.epsilon;
                    }
                }
                // 取平均
                const count = Math.max(1, paramsList.length);
                for (let j = 0; j < numBins; j++) {
                    theoryE[j] /= count;
                    theoryDEdr[j] /= count;
                }
            }

            state.backgroundChartData.energy = {
                centers: radialCenters,
                potential: theoryE,
                dEdr: theoryDEdr,
                localEnergy: {
                    epsilon: avgEpsilon
                }
            };
        } catch (e) {
            console.error('初始化理论能量数据失败:', e);
        }

        console.log('理论曲线数据初始化完成，轨道数:', paramsList.length);

    } catch (error) {
        console.error('初始化理论数据失败:', error);
    }
};

// 后台准备图表数据
window.ElectronCloud.Orbital.updateBackgroundChartData = function () {
    const state = window.ElectronCloud.state;

    if (!window.Hydrogen || state.radialSamples.length === 0) return;

    if (!window.Hydrogen || state.radialSamples.length === 0) return;

    try {
        const angularBins = 180;

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

            // 【参数调整】禁用平滑 (smooth=false)，以便在能量积分中保留采样噪声
            const radialHist = Hydrogen.histogramRadialFromSamples(state.radialSamples, adaptiveBins, dynamicRmax, true, false);

            // 横轴中心
            const centers = new Array(adaptiveBins);
            for (let i = 0; i < adaptiveBins; i++) centers[i] = 0.5 * (radialHist.edges[i] + radialHist.edges[i + 1]);

            // 【关键修复】检测杂化模式，使用正确的物理公式
            const isHybridMode = state.isHybridMode === true;
            const atomType = state.currentAtom || 'H';
            const atomZ = (window.SlaterBasis && window.SlaterBasis[atomType]) ? window.SlaterBasis[atomType].Z : 1;
            const values = new Array(adaptiveBins);
            const wave = new Array(adaptiveBins);

            for (let i = 0; i < adaptiveBins; i++) {
                if (isHybridMode) {
                    // 【杂化模式】使用专门的杂化径向PDF
                    values[i] = Hydrogen.hybridRadialPDF(paramsList, centers[i], atomZ, 1, atomType);
                    let sumWave = 0;
                    const defaultCoeff = 1.0 / Math.sqrt(paramsList.length);
                    for (const params of paramsList) {
                        const coeff = params.coefficient !== undefined ? params.coefficient : defaultCoeff;
                        sumWave += coeff * Math.abs(Hydrogen.radialR(params.n, params.l, centers[i], atomZ, 1, atomType));
                    }
                    wave[i] = sumWave;
                } else {
                    // 【普通模式】各轨道PDF的平均
                    let sumPDF = 0;
                    let sumWave = 0;
                    for (const params of paramsList) {
                        sumPDF += Hydrogen.radialPDF(params.n, params.l, centers[i], atomZ, 1, atomType);
                        sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, centers[i], atomZ, 1, atomType));
                    }
                    values[i] = sumPDF / paramsList.length;
                    wave[i] = sumWave / paramsList.length;
                }

            }

            state.backgroundChartData.radial = {
                hist: radialHist,
                theory: { centers, values, wave }
            };

            // 准备角向数据
            if (state.angularSamples.length > 0) {
                const angularHist = Hydrogen.histogramThetaFromSamples(state.angularSamples, angularBins, true);

                // 【关键修复】计算角向理论值 - 多轨道时加和
                const centers = new Array(angularBins);
                const values = new Array(angularBins);
                for (let i = 0; i < angularBins; i++) {
                    centers[i] = 0.5 * (angularHist.edges[i] + angularHist.edges[i + 1]);

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

            // ==================== 能量数据预计算 (RHF/Koga-DZ) ====================
            // 【关键修复】使用直方图的 bins 作为 r 坐标网格，确保 counts 与 r 对齐
            if (state.backgroundChartData.radial && state.backgroundChartData.radial.hist) {
                const radialHist = state.backgroundChartData.radial.hist;
                const numBins = radialHist.counts.length;

                // 从直方图边界计算中心点，确保与 counts 对齐
                const radialCenters = new Array(numBins);
                for (let i = 0; i < numBins; i++) {
                    radialCenters[i] = 0.5 * (radialHist.edges[i] + radialHist.edges[i + 1]);
                }

                const atomType = state.currentAtom || 'H';
                const Z = window.SlaterBasis && window.SlaterBasis[atomType] ? window.SlaterBasis[atomType].Z : 1;

                // 1. 准备轨道参数列表 (含对应的原子类型)
                let paramsWithAtom = [];
                if (state.isCompareMode && state.compareMode?.slotConfigs) {
                    paramsWithAtom = state.compareMode.slotConfigs
                        .filter(cfg => cfg.orbital)
                        .map(cfg => {
                            const params = window.Hydrogen.orbitalParamsFromKey(cfg.orbital);
                            return params ? { ...params, atomType: cfg.atom || 'H' } : null;
                        })
                        .filter(Boolean);
                } else {
                    paramsWithAtom = paramsList.map(p => ({ ...p, atomType }));
                }

                // 2. 计算理论能量曲线 (简化版：仅使用 ε·P(r))
                const theoryE = new Float32Array(numBins);
                const theoryDEdr = new Float32Array(numBins);
                let avgEpsilon = -0.5;  // 默认本征值

                if (paramsWithAtom.length > 0 && window.Hydrogen.calculateCumulativeOrbitalEnergy) {
                    for (const p of paramsWithAtom) {
                        const pZ = window.SlaterBasis[p.atomType] ? window.SlaterBasis[p.atomType].Z : 1;
                        const res = window.Hydrogen.calculateCumulativeOrbitalEnergy(
                            p.n, p.l, pZ, p.atomType, radialCenters
                        );
                        if (res && res.E) {
                            for (let j = 0; j < numBins && j < res.E.length; j++) {
                                theoryE[j] += res.E[j];
                                theoryDEdr[j] += (res.dEdr ? res.dEdr[j] : 0);
                            }
                            // 使用返回的 epsilon
                            if (res.epsilon !== undefined) {
                                avgEpsilon = res.epsilon;
                            }
                        }
                    }
                    // 取平均
                    const count = paramsWithAtom.length;
                    for (let j = 0; j < numBins; j++) {
                        theoryE[j] /= count;
                        theoryDEdr[j] /= count;
                    }
                }

                // 3. 基于采样数据计算能量直方图（复用径向采样数据 × ε）
                // 能量密度直方图：radialHist.counts × avgEpsilon （采样直方图 × 本征值）
                const samplingDensity = new Float32Array(numBins);
                for (let j = 0; j < numBins; j++) {
                    samplingDensity[j] = radialHist.counts[j] * avgEpsilon;
                }

                // 能量累计直方图：累加 samplingDensity
                const samplingCumulative = new Float32Array(numBins);
                let cumSum = 0;
                const dr = radialHist.dr || (radialCenters[1] - radialCenters[0]);
                for (let j = 0; j < numBins; j++) {
                    cumSum += samplingDensity[j] * dr;
                    samplingCumulative[j] = cumSum;
                }

                state.backgroundChartData.energy = {
                    centers: radialCenters,
                    edges: radialHist.edges,
                    // 采样数据（直方图显示用）
                    samplingDensity: samplingDensity,     // 采样 × ε
                    samplingCumulative: samplingCumulative, // 累计
                    // 理论数据（曲线显示用）
                    theoryDensity: theoryDEdr,            // 理论 ε·P(r)
                    theoryCumulative: theoryE,            // 理论 E(R)
                    epsilon: avgEpsilon
                };

                console.log('能量数据已预计算, 轨道数:', paramsWithAtom.length, ', ε =', avgEpsilon.toFixed(4), 'Ha');
            }

            console.log('后台图表数据已更新，径向数据点:', state.radialSamples.length, '轨道数:', paramsList.length);

            console.log('后台图表数据已更新，径向数据点:', state.radialSamples.length, '轨道数:', paramsList.length);
        } // Close if (state.radialSamples.length > 0)

    } catch (error) {
        console.error('后台更新图表数据失败:', error);
    }
};

// 绘制概率图表
window.ElectronCloud.Orbital.drawProbabilityChart = function (final = true) {
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

    // 确保能量图有可用的纯理论数据（不依赖采样）
    if (['potential', 'dEdr', 'localEnergy'].includes(type) && (!state.backgroundChartData || !state.backgroundChartData.energy)) {
        if (window.ElectronCloud?.Orbital?.initTheoryData) {
            window.ElectronCloud.Orbital.initTheoryData();
        }
    }

    // 检查是否是对比模式
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;

    // 【关键修复】对比模式下图表必须优先走 compare 渲染路径；否则会落入“多选/平均”的后台数据渲染，表现成多选模式
    // 即使采样尚未开始（orbitalSamplesMap 为空），compare 渲染也应能显示纯理论曲线（由 DataPanel.renderChartCompare 兜底）
    if (isCompareMode && ['radial', 'angular', 'phi', 'azimuthal', 'energyDensity', 'energyCumulative'].includes(type)) {
        let chartData = state.orbitalSamplesMap || {};

        // 【性能优化】如果是滚动生成过程中的动态更新（!final），只使用最近的数据切片
        // 避免全量重绘导致卡顿（全量数据可能包含百万个点）
        if (!final && chartData && Object.keys(chartData).length > 0) {
            const sliced = {};
            // 只取每条轨道最近的 20000 个点，足够显示分布趋势且极快
            const SLICE_SIZE = 20000;
            for (const key in chartData) {
                const samples = chartData[key];
                if (samples && samples.length > SLICE_SIZE) {
                    sliced[key] = samples.slice(-SLICE_SIZE);
                } else {
                    sliced[key] = samples;
                }
            }
            chartData = sliced;
        } else if (Object.keys(chartData).length > 0) {
            console.log('使用对比模式散点图渲染，轨道数量:', Object.keys(chartData).length);
        }

        DataPanel.renderChartCompare(chartData, type);
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
    } else if (type === 'phi' && state.backgroundChartData.phi) {
        console.log('使用后台准备的φ角向图表数据进行渲染');
        DataPanel.renderChartPhi(state.backgroundChartData.phi.hist, state.backgroundChartData.phi.theory);
        return;
    } else if (type === 'energyDensity' && state.backgroundChartData.energy) {
        // 能量期望密度：采样数据 × ε 作为直方图，理论曲线叠加
        const e = state.backgroundChartData.energy;
        console.log('使用径向采样数据×ε渲染能量密度图');
        DataPanel.renderChartEnergyDensity(
            { counts: Array.from(e.samplingDensity), edges: Array.from(e.edges) },
            { centers: Array.from(e.centers), values: Array.from(e.theoryDensity) }
        );
        return;
    } else if (type === 'energyCumulative' && state.backgroundChartData.energy) {
        // 能量期望累计：累计采样数据作为直方图，理论曲线叠加
        const e = state.backgroundChartData.energy;
        console.log('使用累计采样数据渲染能量累计图');
        DataPanel.renderChartEnergyCumulative(
            { counts: Array.from(e.samplingCumulative), edges: Array.from(e.edges) },
            { centers: Array.from(e.centers), values: Array.from(e.theoryCumulative) }
        );
        return;
    } else if (type === 'potential' || type === 'dEdr' || type === 'localEnergy') {
        // 能量相关图表：如果没有在上面处理（说明没有backgroundData），也不应该 fall through 到重新计算逻辑
        // 重新计算逻辑暂时只支持 radial/angular/phi。
        // 事实上，backgroundChartData.energy 应该总是有的，除非初始化失败。
        if (state.backgroundChartData.energy) {
            const e = state.backgroundChartData.energy;
            // 注意：potential 等图表主要依靠 renderChartPotential 等函数，它们只需要 theory 数据
            // 这里的 dispatch 需要补全
            if (type === 'potential') {
                DataPanel.renderChartPotential({ centers: Array.from(e.centers), values: Array.from(e.potential) });
            } else if (type === 'dEdr') {
                DataPanel.renderChartDEdr({ centers: Array.from(e.centers), values: Array.from(e.dEdr) });
            } else if (type === 'localEnergy') {
                DataPanel.renderChartLocalEnergy({ centers: Array.from(e.centers), values: new Array(e.centers.length).fill(e.localEnergy.epsilon) });
            }
            return;
        }
    }

    // 如果没有后台数据，则重新计算
    console.log('后台数据不可用，重新计算图表数据');
    const angularBins = 180;

    // 【关键修复】使用所有选中的轨道计算理论曲线，而不是单个轨道
    // 在多选模式下，如果未选中任何轨道，不应退化为 '1s'，而应显示空图表
    const isMultiselectMode = ui.multiselectToggle && ui.multiselectToggle.checked;
    // isCompareMode 已在上方定义 (line 936)

    let orbitals = state.currentOrbitals || [];

    // 如果列表为空且不是多选/比照模式（即单选模式），由于必须有一个选中，通常 state.currentOrbital 会有值
    // 但如果用户确实清空了选择（在多选/比照模式允许），我们应该尊重"空选择"
    if (orbitals.length === 0) {
        if (!isMultiselectMode && !isCompareMode) {
            // 单选模式Fallback
            orbitals = [state.currentOrbital || '1s'];
        } else {
            // 多选/比照模式允许空选择：保持空数组
            orbitals = [];
        }
    }

    const paramsList = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    // 如果没有选中轨道，且是多选/比照模式，我们需要渲染一个"空图表"来清除之前的显示
    // 而不是直接 return，否则旧图表会作为"幽灵"保留在屏幕上
    if (paramsList.length === 0) {
        if (type === 'angular') {
            const hist = Hydrogen.histogramThetaFromSamples([], angularBins, true);
            DataPanel.renderChartAngular(hist, null);
        } else if (type === 'phi' || type === 'azimuthal') {
            const hist = Hydrogen.histogramPhiFromSamples([], 180, true);
            DataPanel.renderChartPhi(hist, null);
        }
        // 径向和能量图表可以保留 return，或者也清除。用户反馈主要是角向 ghost line。
        // 为了一致性，我们也清除径向
        if (type === 'radial') {
            const hist = Hydrogen.histogramRadialFromSamples([], 240, 10, true);
            DataPanel.renderChartRadial(hist, null);
        }
        return;
    }

    if (type === 'radial' || type === 'energyDensity' || type === 'energyCumulative') {
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

        // 【关键修复】检测杂化模式，使用正确的物理公式
        const isHybridMode = state.isHybridMode === true;
        const atomType = state.currentAtom || 'H';
        const atomZ = (window.SlaterBasis && window.SlaterBasis[atomType]) ? window.SlaterBasis[atomType].Z : 1;
        const values = new Array(adaptiveBins);
        const wave = new Array(adaptiveBins);

        for (let i = 0; i < adaptiveBins; i++) {
            if (isHybridMode) {
                // 【杂化模式】使用专门的杂化径向PDF
                values[i] = Hydrogen.hybridRadialPDF(paramsList, centers[i], atomZ, 1, atomType);
                let sumWave = 0;
                const defaultCoeff = 1.0 / Math.sqrt(paramsList.length);
                for (const params of paramsList) {
                    const coeff = params.coefficient !== undefined ? params.coefficient : defaultCoeff;
                    sumWave += coeff * Math.abs(Hydrogen.radialR(params.n, params.l, centers[i], atomZ, 1, atomType));
                }
                wave[i] = sumWave;
            } else {
                // 【普通模式】各轨道PDF的平均
                let sumPDF = 0;
                let sumWave = 0;
                for (const params of paramsList) {
                    sumPDF += Hydrogen.radialPDF(params.n, params.l, centers[i], atomZ, 1, atomType);
                    sumWave += Math.abs(Hydrogen.radialR(params.n, params.l, centers[i], atomZ, 1, atomType));
                }
                values[i] = sumPDF / paramsList.length;
                wave[i] = sumWave / paramsList.length;
            }
        }

        // 绘制直方图和理论曲线
        if (type === 'radial') {
            DataPanel.renderChartRadial(hist, { centers, values, wave });
        }
        // 能量图表不再使用 fallback - 必须依赖预计算的 backgroundChartData.energy
    } else if (type === 'angular') {
        // θ角向分布
        const hist = Hydrogen.histogramThetaFromSamples(state.angularSamples, angularBins, true);

        // 【关键修复】计算角向理论值 - 多轨道时加和
        const centers = new Array(angularBins);
        const values = new Array(angularBins);
        for (let i = 0; i < angularBins; i++) {
            centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);

            let sumAngularPDF = 0;
            for (const params of paramsList) {
                sumAngularPDF += Hydrogen.angularPDF_Theta(params.angKey.l, params.angKey.m, centers[i]);
            }
            values[i] = sumAngularPDF / paramsList.length;
        }

        DataPanel.renderChartAngular(hist, { centers, values });
    } else if (type === 'phi' || type === 'azimuthal') {
        // φ角向分布 (azimuthal)
        const phiBins = 180;
        const hist = Hydrogen.histogramPhiFromSamples(state.phiSamples, phiBins, true);

        // 【修复】φ方向理论值应根据轨道的m值和类型（cos/sin）计算
        // 对于m=0：均匀分布(1/2π)
        // 对于m≠0 cos型：cos²(mφ)/π
        // 对于m≠0 sin型：sin²(mφ)/π
        const centers = new Array(phiBins);
        const values = new Array(phiBins);
        for (let i = 0; i < phiBins; i++) {
            centers[i] = 0.5 * (hist.edges[i] + hist.edges[i + 1]);

            // 多轨道时取平均
            let sumPhiPDF = 0;
            for (const params of paramsList) {
                sumPhiPDF += Hydrogen.angularPDF_Phi(params.angKey.m, params.angKey.t, centers[i]);
            }
            values[i] = sumPhiPDF / paramsList.length;
        }

        DataPanel.renderChartPhi(hist, { centers, values });
    } else {
        console.warn('未知的图表类型或无法动态计算:', type);
    }
};

/**
 * 更新比照模式下某个slot的可见性
 * 【重构】使用slot索引控制可见性，而非轨道键
 * @param {number} orbitalIndex - 在activeSlots数组中的索引（Worker中的orbitalIndex）
 * @param {boolean} visible - 是否可见
 */
window.ElectronCloud.Orbital.updateCompareSlotVisibility = function (orbitalIndex, visible) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    console.log(`updateCompareSlotVisibility: orbitalIndex = ${orbitalIndex}, visible = ${visible} `);
    console.log(`state.pointCount = ${state.pointCount}, state.renderingCompleted = ${state.renderingCompleted} `);
    console.log(`pointOrbitalIndices存在？${!!state.pointOrbitalIndices} `);

    // 只有渲染完成后才能切换可见性
    if (!state.points || !state.renderingCompleted || !ui.compareToggle.checked) {
        console.log('渲染未完成或不在比照模式，返回');
        return;
    }

    const positions = state.points.geometry.attributes.position.array;

    // 如果还没有备份原始位置，先备份
    if (!state.originalPositions) {
        state.originalPositions = new Float32Array(positions);
        console.log('已备份原始位置');
    }

    // 统计匹配的点数
    let matchCount = 0;

    // 遍历所有点，找到属于这个slot的点并切换可见性
    for (let i = 0; i < state.pointCount; i++) {
        const pointOrbitalIndex = state.pointOrbitalIndices ? state.pointOrbitalIndices[i] : 0;

        // 检查这个点是否属于当前轨道索引
        if (pointOrbitalIndex === orbitalIndex) {
            matchCount++;
            const posIdx = i * 3;
            if (visible) {
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
    }

    state.points.geometry.attributes.position.needsUpdate = true;
    console.log(`orbitalIndex = ${orbitalIndex} 可见性已更新为: ${visible}, 匹配${matchCount} 个点`);

    // 【关键新增】同步更新对应的等值面轮廓可见性
    if (window.ElectronCloud.Visualization && window.ElectronCloud.Visualization.updateCompareContourVisibility) {
        window.ElectronCloud.Visualization.updateCompareContourVisibility(orbitalIndex, visible);
    }
};
