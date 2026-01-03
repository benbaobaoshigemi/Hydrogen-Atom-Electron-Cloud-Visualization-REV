// 点采样逻辑模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Sampling = {};

// ==================== Web Worker 管理 ====================

// Worker 池配置
const WORKER_POOL_SIZE = navigator.hardwareConcurrency || 4; // 使用 CPU 核心数
const workerPool = [];
let workerReady = 0;
let taskIdCounter = 0;
let useWorkers = true; // 是否使用 Worker（可降级）
let workersInitialized = false; // 防止重复初始化
let lastSubmitTime = 0; // 上次提交时间（节流）
let samplingSessionId = 0; // 采样会话 ID，用于防止旧 Worker 结果污染新会话

// 初始化 Worker 池
window.ElectronCloud.Sampling.initWorkers = function () {
    // 防止重复初始化
    if (workersInitialized) {
        return;
    }
    workersInitialized = true;

    if (typeof Worker === 'undefined') {
        console.warn('Web Workers 不可用，将使用主线程采样');
        useWorkers = false;
        return;
    }

    try {
        for (let i = 0; i < WORKER_POOL_SIZE; i++) {
            const worker = new Worker('sampling-worker.js');

            worker.onmessage = function (e) {
                const { type, taskId, sessionId, result } = e.data;

                if (type === 'ready') {
                    workerReady++;
                    console.log(`Worker ${i} 就绪 (${workerReady}/${WORKER_POOL_SIZE})`);
                    worker.busy = false;
                } else if (type === 'sample-result') {
                    // 处理采样结果
                    window.ElectronCloud.Sampling.handleWorkerResult(result, taskId, sessionId);
                    worker.busy = false;
                }
            };

            worker.onerror = function (err) {
                console.error('Worker 错误:', err);
                worker.busy = false;
            };

            worker.busy = false;
            workerPool.push(worker);
        }

        console.log(`初始化 ${WORKER_POOL_SIZE} 个 Web Worker 用于并行采样`);
    } catch (err) {
        console.warn('Worker 初始化失败，将使用主线程采样:', err);
        useWorkers = false;
    }
};

// 处理 Worker 返回的结果
window.ElectronCloud.Sampling.handleWorkerResult = function (result, taskId, sessionId) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    // 检查 session ID，忽略来自旧会话的结果
    // 注意：sessionId 可能是 undefined（旧 Worker 或首次初始化），此时允许通过
    if (sessionId !== undefined && sessionId !== samplingSessionId) {
        console.log(`忽略旧会话的 Worker 结果 (session ${sessionId} vs current ${samplingSessionId})`);
        return;
    }

    if (!state.points || !result.points.length) return;

    // 【关键修复】确保采样索引映射存在：用于让 orbitalSamplesMap 与点云一一对应
    if (!state.orbitalSampleIndexByPoint) {
        state.orbitalSampleIndexByPoint = {};
    }

    const positions = state.points.geometry.attributes.position.array;
    const colorsAttr = state.points.geometry.getAttribute('color');
    const colors = colorsAttr ? colorsAttr.array : null;

    // 确保 pointOrbitalIndices 数组存在
    if (!state.pointOrbitalIndices || state.pointOrbitalIndices.length < state.MAX_POINTS) {
        state.pointOrbitalIndices = new Int16Array(state.MAX_POINTS);
    }

    // 将 Worker 采样的点添加到主线程
    const startPointIndex = state.pointCount;
    const maxAdd = Math.max(0, state.MAX_POINTS - startPointIndex);
    const pointsToAdd = Math.min(result.points.length, maxAdd);

    for (let pIdx = 0; pIdx < pointsToAdd; pIdx++) {
        const point = result.points[pIdx];
        const pointIndex = state.pointCount;

        const i3 = pointIndex * 3;
        positions[i3] = point.x;
        positions[i3 + 1] = point.y;
        positions[i3 + 2] = point.z;

        if (colors) {
            colors[i3] = point.r;
            colors[i3 + 1] = point.g;
            colors[i3 + 2] = point.b;
        }

        // 同步 baseColors
        if (state.baseColors) {
            state.baseColors[i3] = point.r;
            state.baseColors[i3 + 1] = point.g;
            state.baseColors[i3 + 2] = point.b;
            state.baseColorsCount = state.pointCount + 1;
        }

        // 存储每个点的轨道索引（用于多轨道模式的等值面计算）
        state.pointOrbitalIndices[pointIndex] = point.orbitalIndex >= 0 ? point.orbitalIndex : 0;

        // 记录轨道点映射（用于比照模式开关）
        if (point.orbitalIndex >= 0) {
            // 【关键修复】比照模式下 orbitalPointsMap 也必须使用唯一键（原子+轨道+slot）
            let orbitalKey = state.currentOrbitals[point.orbitalIndex];
            const isCompareModeForPoint = ui.compareToggle && ui.compareToggle.checked;
            if (isCompareModeForPoint && state.compareMode && state.compareMode.activeSlots) {
                const slotConfig = state.compareMode.activeSlots[point.orbitalIndex];
                if (slotConfig) {
                    orbitalKey = `${slotConfig.atom || 'H'}_${slotConfig.orbital}_slot${slotConfig.slotIndex}`;
                }
            }

            if (orbitalKey) {
                if (!state.orbitalPointsMap[orbitalKey]) {
                    state.orbitalPointsMap[orbitalKey] = [];
                }
                state.orbitalPointsMap[orbitalKey].push(pointIndex);
            }
        }

        state.pointCount++;
    }

    // 添加采样数据（保证与 pointIndex 对齐）
    const samplesToAdd = Math.min(result.samples.length, pointsToAdd);
    for (let sIdx = 0; sIdx < samplesToAdd; sIdx++) {
        const sample = result.samples[sIdx];
        const pointIndex = startPointIndex + sIdx;

        state.radialSamples.push(sample.r);
        state.angularSamples.push(sample.theta);
        state.phiSamples.push(sample.phi);

        if (sample.orbitalKey) {
            // 【关键修复】比照模式下使用slot索引+原子类型作为唯一键，避免同轨道覆盖
            const ui = window.ElectronCloud.ui;
            const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
            let sampleKey = sample.orbitalKey;

            if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                // 使用orbitalIndex找到对应的slot配置
                const slotConfig = state.compareMode.activeSlots[sample.orbitalIndex];
                console.log(`比照模式数据键构建: orbitalIndex=${sample.orbitalIndex}, activeSlots长度=${state.compareMode.activeSlots.length}, slotConfig=`, slotConfig);
                if (slotConfig) {
                    // 使用"原子_轨道_slot索引"作为唯一键
                    sampleKey = `${slotConfig.atom || 'H'}_${sample.orbitalKey}_slot${slotConfig.slotIndex}`;
                    console.log(`sampleKey=${sampleKey}`);
                }
            }

            if (!state.orbitalSamplesMap[sampleKey]) {
                state.orbitalSamplesMap[sampleKey] = [];
            }
            const arr = state.orbitalSamplesMap[sampleKey];
            arr.push({
                r: sample.r,
                theta: sample.theta,
                phi: sample.phi,
                probability: 0, // Worker 中已计算过
                orbitalIndex: sample.orbitalIndex, // 保存slot索引用于后续处理
                pointIndex: pointIndex
            });

            // 记录 pointIndex -> (key, idx) 映射，便于滚动更新时精确删除对应样本
            state.orbitalSampleIndexByPoint[pointIndex] = { key: sampleKey, idx: arr.length - 1 };
        }
    }

    // 更新最远距离
    if (result.farthestDistance > state.farthestDistance) {
        state.farthestDistance = result.farthestDistance;
        state.samplingBoundary = Math.max(state.samplingBoundary, state.farthestDistance * 1.05);
    }

    // 更新几何体
    state.points.geometry.setDrawRange(0, state.pointCount);
    state.points.geometry.attributes.position.needsUpdate = true;
    if (colorsAttr) colorsAttr.needsUpdate = true;

    // 更新计数显示
    const ui_pointCountSpan = window.ElectronCloud.ui.pointCountSpan;
    if (ui_pointCountSpan) {
        ui_pointCountSpan.textContent = state.pointCount;
    }
};

// 提交采样任务到 Worker 池
window.ElectronCloud.Sampling.submitSamplingTask = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const constants = window.ElectronCloud.constants;

    if (state.pointCount >= state.MAX_POINTS) return;

    // 检查是否有轨道可采样
    if (!state.currentOrbitals || state.currentOrbitals.length === 0) return;

    // 检查是否有空闲 Worker
    const availableWorkers = workerPool.filter(w => !w.busy);
    if (availableWorkers.length === 0) return;

    // 【节流机制】根据速度滑条控制提交频率
    // 滑条范围 20000-800000，值越大，提交越频繁
    const now = performance.now();
    const speedValue = window.ElectronCloud.Sampling.getAttemptsPerFrame();
    // 将速度值映射到提交间隔：低速 = 长间隔，高速 = 短间隔
    // 速度 20000 → 间隔 50ms，速度 800000 → 间隔 0ms（无限制）
    const minInterval = Math.max(0, 50 * (1 - (speedValue - 20000) / (800000 - 20000)));

    if (now - lastSubmitTime < minInterval) {
        return; // 节流：跳过此次提交
    }
    lastSubmitTime = now;

    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    const isMultiselectMode = ui.multiselectToggle && ui.multiselectToggle.checked;
    const isHybridMode = state.isHybridMode; // 【新增】杂化模式标志
    const phaseOn = state.usePhaseColoring || false;

    // 准备比照颜色（只传递值，不传递对象）
    // 【关键修复】比照模式颜色必须绑定到 UI slotIndex，避免空 slot 导致颜色错位
    let compareColors = [];
    if (constants.compareColors) {
        if (isCompareMode && state.compareMode && Array.isArray(state.compareMode.activeSlots) && state.compareMode.activeSlots.length > 0) {
            const palette = constants.compareColors;
            compareColors = state.compareMode.activeSlots.map(slot => {
                const idx = (slot && Number.isInteger(slot.slotIndex)) ? slot.slotIndex : 0;
                const safe = ((idx % palette.length) + palette.length) % palette.length;
                return palette[safe].value;
            });
        } else {
            compareColors = constants.compareColors.map(c => c.value);
        }
    }

    const taskId = ++taskIdCounter;
    // 根据速度滑条调整每个 Worker 的工作量
    const totalAttempts = speedValue;
    const attemptsPerWorker = Math.ceil(totalAttempts / availableWorkers.length);
    const pointsPerWorker = Math.ceil(state.pointsPerFrame / availableWorkers.length);

    // 为每个空闲 Worker 分配任务
    for (let i = 0; i < availableWorkers.length; i++) {
        const worker = availableWorkers[i];
        const task = {
            type: 'sample',
            taskId: taskId + i,
            sessionId: samplingSessionId, // 添加会话 ID
            data: {
                orbitals: state.currentOrbitals,
                samplingBoundary: state.samplingBoundary,
                maxAttempts: attemptsPerWorker,
                targetPoints: pointsPerWorker,
                isIndependentMode: !isHybridMode && (isCompareMode || isMultiselectMode),
                isMultiselectMode: !isHybridMode && isMultiselectMode,
                isCompareMode: isCompareMode, // 【新增】明确传递比照模式标志
                isHybridMode: isHybridMode,
                phaseOn: phaseOn,
                compareColors: compareColors,
                atomType: state.currentAtom || 'H',
                // 【重构】比照模式使用slotConfigs，每个slot独立处理
                slotConfigs: isCompareMode && state.compareMode ? state.compareMode.slotConfigs : null
            }
        };

        worker.busy = true;
        worker.postMessage(task);
    }
};

// 重置采样会话（使旧的 Worker 结果无效）
window.ElectronCloud.Sampling.resetSamplingSession = function () {
    samplingSessionId++;
    console.log(`采样会话重置为 ${samplingSessionId}`);
};

// 检查是否使用 Worker
window.ElectronCloud.Sampling.isUsingWorkers = function () {
    // 【关键修复】杂化模式下禁用 Worker
    // 因为 Worker 无法访问主线程的 hybridRenderMode 状态
    // 且 Worker 中的杂化采样函数不支持"全部"模式
    const state = window.ElectronCloud.state;
    if (state && state.isHybridMode) {
        return false; // 杂化模式强制使用主线程采样
    }
    return useWorkers && workerReady > 0;
};

// ==================== 原有采样逻辑（作为降级方案）====================

// 是否启用重要性采样（性能优化模式）
let useImportanceSampling = true;

// 更新点的位置（主要采样逻辑 - 支持重要性采样）
window.ElectronCloud.Sampling.updatePoints = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const constants = window.ElectronCloud.constants;

    if (!state.points) return;

    // 【关键修复】检查是否使用 Worker
    // 注意：isUsingWorkers() 在杂化模式下会强制返回 false
    if (window.ElectronCloud.Sampling.isUsingWorkers()) {
        window.ElectronCloud.Sampling.submitSamplingTask();
        return;
    }

    const positions = state.points.geometry.attributes.position.array;
    const colorsAttr = state.points.geometry.getAttribute('color') || (function () {
        const colors = new Float32Array(state.points.geometry.attributes.position.count * 3);
        state.points.geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
        return state.points.geometry.getAttribute('color');
    })();

    // 确保 pointOrbitalIndices 数组存在（用于多轨道模式的等值面计算）
    if (!state.pointOrbitalIndices || state.pointOrbitalIndices.length < state.MAX_POINTS) {
        state.pointOrbitalIndices = new Int16Array(state.MAX_POINTS);
    }

    const paramsList = state.currentOrbitals.map(k => Hydrogen.orbitalParamsFromKey(k)).filter(Boolean);
    if (!paramsList.length) return;

    let newPoints = 0;
    const attemptPerFrame = window.ElectronCloud.Sampling.getAttemptsPerFrame();
    const samplingVolumeSize = state.samplingBoundary * 2;

    // 性能优化：允许更长的计算时间，充分利用帧间隔
    const startTime = performance.now();
    const maxTimePerFrame = 15; // 毫秒（60fps下每帧约16.67ms）

    // 判断是否使用独立模式（多选/比照）- 但杂化模式除外
    const isHybridMode = state.isHybridMode;
    const isIndependentMode = !isHybridMode && ((ui.compareToggle && ui.compareToggle.checked) ||
        (ui.multiselectToggle && ui.multiselectToggle.checked));
    const isMultiselectMode = !isHybridMode && ui.multiselectToggle && ui.multiselectToggle.checked;

    let attempts = 0;
    while (attempts < attemptPerFrame && newPoints < state.pointsPerFrame && state.pointCount < state.MAX_POINTS) {
        // 每1000次尝试检查一次时间（减少开销）
        if (attempts % 1000 === 0 && attempts > 0 && (performance.now() - startTime) > maxTimePerFrame) {
            break;
        }
        attempts++;

        // 【杂化模式】使用专门的杂化采样逻辑
        if (isHybridMode) {
            const success = window.ElectronCloud.Sampling.processHybridModePoint(
                paramsList, positions, colorsAttr.array
            );
            if (success) {
                newPoints++;
            }
            continue;
        }

        // 根据采样模式选择不同的策略
        if (useImportanceSampling && window.Hydrogen.importanceSample) {
            // 重要性采样模式
            const success = window.ElectronCloud.Sampling.processImportanceSamplingPoint(
                paramsList, positions, colorsAttr.array, isIndependentMode, isMultiselectMode
            );
            if (success) {
                newPoints++;
            }
        } else {
            // 传统拒绝采样模式（降级方案）
            const x = (Math.random() - 0.5) * samplingVolumeSize;
            const y = (Math.random() - 0.5) * samplingVolumeSize;
            const z = (Math.random() - 0.5) * samplingVolumeSize;
            const r = Math.sqrt(x * x + y * y + z * z);

            if (r === 0) continue;

            const theta = Math.acos(z / r);
            const phi = Math.atan2(y, x);

            if (isIndependentMode) {
                const success = window.ElectronCloud.Sampling.processIndependentModePoint(
                    x, y, z, r, theta, phi, paramsList, positions, colorsAttr.array, isMultiselectMode
                );
                if (success) {
                    newPoints++;
                }
            } else {
                const success = window.ElectronCloud.Sampling.processNormalModePoint(
                    x, y, z, r, theta, phi, paramsList, positions, colorsAttr.array
                );
                if (success) {
                    newPoints++;
                }
            }
        }
    }

    // 采样完成时使用实际数据更新角向分布
    if (state.samplingCompleted && ui.angular3dToggle && ui.angular3dToggle.checked) {
        if (!state.angularUpdated) {
            state.angularUpdated = true;
            setTimeout(() => {
                window.ElectronCloud.Visualization.updateAngularOverlayFromSamples();
                state.angularUpdated = false;
            }, 100);
        }
    } else if (ui.angular3dToggle && ui.angular3dToggle.checked && state.pointCount % 10000 === 0) {
        window.ElectronCloud.Visualization.updateAngularOverlayFromSamples();
    }

    state.points.geometry.setDrawRange(0, state.pointCount);
    state.points.geometry.attributes.position.needsUpdate = true;
    if (state.points.geometry.attributes.color) state.points.geometry.attributes.color.needsUpdate = true;

    const ui_pointCountSpan = window.ElectronCloud.ui.pointCountSpan;
    if (ui_pointCountSpan) {
        ui_pointCountSpan.textContent = state.pointCount;
    }
};

// 处理独立模式下的点采样（多选模式和比照模式共用）
// 【重要修复】使用轮流采样（Round-Robin）替代随机选择
// isMultiselectMode: 是否为多选模式，用于决定颜色方案
window.ElectronCloud.Sampling.processIndependentModePoint = function (x, y, z, r, theta, phi, paramsList, positions, colors, isMultiselectMode) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    const ui = window.ElectronCloud.ui;

    // 【修复】使用轮流选择而不是随机选择，确保每个轨道获得公平的采样机会
    const orbitalIndex = state.roundRobinIndex % paramsList.length;
    state.roundRobinIndex++;
    const p = paramsList[orbitalIndex];

    // 【关键修复】获取当前轨道对应的原子类型（比照模式支持不同原子）
    const orbitalKey = state.currentOrbitals[orbitalIndex];
    let currentAtom = state.currentAtom || 'H';
    // 【修复】优先使用基于索引的原子列表（支持同名轨道对比不同原子），降级使用 Map
    if (state.compareMode && state.compareMode.atoms && state.compareMode.atoms[orbitalIndex]) {
        currentAtom = state.compareMode.atoms[orbitalIndex];
    } else if (state.compareMode && state.compareMode.orbitalAtomMap && state.compareMode.orbitalAtomMap[orbitalKey]) {
        currentAtom = state.compareMode.orbitalAtomMap[orbitalKey];
    }
    const density = Hydrogen.density3D_real(p.angKey, p.n, p.l, r, theta, phi, 1, 1, currentAtom);

    // 【关键】为每个轨道独立计算scaleFactor
    // 确保每个轨道按照自己的概率密度分布采样，而不是混在一起
    // 这样可以保证各轨道的边界完全独立
    const scaleFactor = Math.min(200, 3.0 * Math.pow(p.n, 3));

    if (Math.random() < density * scaleFactor) {
        const i3 = state.pointCount * 3;
        positions[i3] = x;
        positions[i3 + 1] = y;
        positions[i3 + 2] = z;

        let r_color, g_color, b_color;

        if (isMultiselectMode) {
            // 多选模式：支持相位显示功能
            const phaseOn = window.ElectronCloud.state.usePhaseColoring;
            let sign = 0;
            if (phaseOn) {
                // 只计算当前选中轨道的波函数，不叠加其他轨道
                const R = Hydrogen.radialR(p.n, p.l, r, 1, 1, currentAtom);
                const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                const psi = R * Y;
                sign = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
            }

            // 使用全局 phaseColors 常量
            const pc = window.ElectronCloud.constants.phaseColors;
            if (sign > 0) {
                r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
            } else if (sign < 0) {
                r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
            } else {
                r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
            }
        } else {
            // 比照模式：使用固定颜色区分不同轨道
            // 【关键修复】颜色绑定到 UI slotIndex，避免空 slot 造成颜色前移
            let colorIndex = orbitalIndex;
            const ui = window.ElectronCloud.ui;
            const isCompareModeNow = ui && ui.compareToggle && ui.compareToggle.checked;
            if (isCompareModeNow && state.compareMode && Array.isArray(state.compareMode.activeSlots)) {
                const slot = state.compareMode.activeSlots[orbitalIndex];
                if (slot && Number.isInteger(slot.slotIndex)) colorIndex = slot.slotIndex;
            }

            if (constants.compareColors && constants.compareColors.length > 0) {
                const safe = ((colorIndex % constants.compareColors.length) + constants.compareColors.length) % constants.compareColors.length;
                const color = constants.compareColors[safe].value;
                r_color = color[0];
                g_color = color[1];
                b_color = color[2];
            } else {
                r_color = 1;
                g_color = 1;
                b_color = 1;
            }
        }

        colors[i3] = r_color;
        colors[i3 + 1] = g_color;
        colors[i3 + 2] = b_color;

        // 同步更新baseColors（用于星空闪烁模式）
        if (state.baseColors) {
            state.baseColors[i3] = r_color;
            state.baseColors[i3 + 1] = g_color;
            state.baseColors[i3 + 2] = b_color;
            state.baseColorsCount = state.pointCount + 1;
        }

        // 记录这个点属于哪个轨道（用于后续的显示开关）
        const orbitalKey = state.currentOrbitals[orbitalIndex];
        if (!state.orbitalPointsMap[orbitalKey]) {
            state.orbitalPointsMap[orbitalKey] = [];
        }
        state.orbitalPointsMap[orbitalKey].push(state.pointCount);

        // 存储轨道索引（用于多轨道模式的等值面计算）
        if (state.pointOrbitalIndices) {
            state.pointOrbitalIndices[state.pointCount] = orbitalIndex;
        }

        // 记录采样数据（用于图表绘制）
        if (!state.orbitalSamplesMap[orbitalKey]) {
            state.orbitalSamplesMap[orbitalKey] = [];
        }
        state.orbitalSamplesMap[orbitalKey].push({
            r: r,
            theta: theta,
            probability: density
        });

        state.pointCount++;
        state.radialSamples.push(r);
        state.angularSamples.push(theta);
        state.phiSamples.push(phi);

        window.ElectronCloud.Sampling.updateFarthestDistance(r);

        return true;
    }

    return false;
};

// 处理正常模式下的点采样
window.ElectronCloud.Sampling.processNormalModePoint = function (x, y, z, r, theta, phi, paramsList, positions, colors) {
    const state = window.ElectronCloud.state;

    let probability = 0;

    // 计算所有轨道的概率密度叠加，同时找出最小n值
    let minN = paramsList[0].n;
    // 获取当前原子类型
    const currentAtom = state.currentAtom || 'H';
    for (let i = 0; i < paramsList.length; i++) {
        const p = paramsList[i];
        const density = Hydrogen.density3D_real(p.angKey, p.n, p.l, r, theta, phi, 1, 1, currentAtom);
        probability += density;
        if (p.n < minN) minN = p.n;
    }

    // 使用最小n值计算scaleFactor，确保低n轨道（高密度区域）不被截断
    // 除以sqrt(轨道数)是因为叠加的密度峰值通常小于各轨道峰值之和
    const scaleFactor = Math.min(200, 3.0 * Math.pow(minN, 3) / Math.sqrt(paramsList.length));

    if (Math.random() < probability * scaleFactor) {
        const i3 = state.pointCount * 3;
        positions[i3] = x;
        positions[i3 + 1] = y;
        positions[i3 + 2] = z;

        // 多选模式和默认模式：都支持相位显示功能
        const phaseOn = window.ElectronCloud.state.usePhaseColoring;
        let sign = 0;
        if (phaseOn) {
            let psi = 0;
            for (const p of paramsList) {
                const R = Hydrogen.radialR(p.n, p.l, r, 1, 1, currentAtom);
                const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                psi += R * Y;
            }
            sign = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
        }

        let r_color, g_color, b_color;
        // 使用全局 phaseColors 常量
        const pc = window.ElectronCloud.constants.phaseColors;
        if (sign > 0) {
            r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
        } else if (sign < 0) {
            r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
        } else {
            r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
        }

        colors[i3] = r_color;
        colors[i3 + 1] = g_color;
        colors[i3 + 2] = b_color;

        // 同步更新baseColors（用于星空闪烁模式）
        if (state.baseColors) {
            state.baseColors[i3] = r_color;
            state.baseColors[i3 + 1] = g_color;
            state.baseColors[i3 + 2] = b_color;
            // 更新计数器（pointCount还没有增加，所以是 pointCount + 1）
            state.baseColorsCount = state.pointCount + 1;
        }

        // 存储轨道索引（正常模式为 0，表示叠加模式/单轨道）
        if (state.pointOrbitalIndices) {
            state.pointOrbitalIndices[state.pointCount] = 0;
        }

        state.pointCount++;
        state.radialSamples.push(r);
        state.angularSamples.push(theta);
        state.phiSamples.push(phi);

        window.ElectronCloud.Sampling.updateFarthestDistance(r);

        return true;
    }

    return false;
};

// 更新最远距离和相关参数
window.ElectronCloud.Sampling.updateFarthestDistance = function (r) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const constants = window.ElectronCloud.constants;

    if (r > state.farthestDistance) {
        state.farthestDistance = r;
        state.samplingBoundary = Math.max(state.samplingBoundary, state.farthestDistance * 1.05);

        // 注意：坐标系的实时更新现在主要在动画循环中处理
        // 这里不再需要显式调用 onAxesSizeChange，因为动画循环会自动处理

        if (state.angularOverlay && ui.angular3dToggle && ui.angular3dToggle.checked) {
            if (state.pointCount % 5000 === 0) {
                window.ElectronCloud.Visualization.updateAngularOverlay();
            }
        }
    }
};

// 获取每帧尝试次数
window.ElectronCloud.Sampling.getAttemptsPerFrame = function () {
    const ui = window.ElectronCloud.ui;

    if (!ui.speedRange) return 20000;

    const v = parseInt(ui.speedRange.value, 10);
    if (isNaN(v) || v <= 0) return 20000;
    return v;
};

// ==================== 重要性采样处理函数 ====================

/**
 * 使用重要性采样方法处理单个采样点
 * 
 * 【重要修复】多选/比照模式使用轮流采样（Round-Robin）策略：
 * - 每个轨道轮流获得采样机会
 * - 每个轨道使用自己的重要性采样提议分布
 * - 避免了随机选择轨道导致的密度偏差
 * 
 * @param {Array} paramsList - 轨道参数列表
 * @param {Float32Array} positions - 位置数组
 * @param {Float32Array} colors - 颜色数组
 * @param {boolean} isIndependentMode - 是否为独立模式（多选/比照）
 * @param {boolean} isMultiselectMode - 是否为多选模式
 * @returns {boolean} - 是否成功添加了一个点
 */
window.ElectronCloud.Sampling.processImportanceSamplingPoint = function (
    paramsList, positions, colors, isIndependentMode, isMultiselectMode
) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;

    // 【修复】使用轮流选择而不是随机选择，确保每个轨道获得公平的采样机会
    let orbitalIndex = 0;
    if (isIndependentMode && paramsList.length > 1) {
        orbitalIndex = state.roundRobinIndex % paramsList.length;
        state.roundRobinIndex++;
    }

    const p = paramsList[orbitalIndex];

    // 【关键修复】获取当前轨道对应的原子类型（比照模式支持不同原子）
    const orbitalKey = state.currentOrbitals[orbitalIndex];
    let currentAtom = state.currentAtom || 'H';
    if (state.compareMode && state.compareMode.orbitalAtomMap && state.compareMode.orbitalAtomMap[orbitalKey]) {
        currentAtom = state.compareMode.orbitalAtomMap[orbitalKey];
    }

    // 使用重要性采样生成点（针对当前轨道优化）
    const result = Hydrogen.importanceSample(p.n, p.l, p.angKey, state.samplingBoundary, currentAtom);

    if (!result || !result.accepted) {
        return false;
    }

    const { x, y, z, r, theta, phi } = result;

    // 边界检查
    if (r > state.samplingBoundary * 2) {
        return false;
    }

    const i3 = state.pointCount * 3;
    positions[i3] = x;
    positions[i3 + 1] = y;
    positions[i3 + 2] = z;

    // 计算颜色
    let r_color, g_color, b_color;

    if (isIndependentMode && !isMultiselectMode) {
        // 比照模式：使用固定颜色区分不同轨道
        // 【关键修复】颜色绑定到 UI slotIndex，避免空 slot 造成颜色前移
        let colorIndex = orbitalIndex;
        const ui = window.ElectronCloud.ui;
        const isCompareModeNow = ui && ui.compareToggle && ui.compareToggle.checked;
        if (isCompareModeNow && state.compareMode && Array.isArray(state.compareMode.activeSlots)) {
            const slot = state.compareMode.activeSlots[orbitalIndex];
            if (slot && Number.isInteger(slot.slotIndex)) colorIndex = slot.slotIndex;
        }

        if (constants.compareColors && constants.compareColors.length > 0) {
            const safe = ((colorIndex % constants.compareColors.length) + constants.compareColors.length) % constants.compareColors.length;
            const color = constants.compareColors[safe].value;
            r_color = color[0];
            g_color = color[1];
            b_color = color[2];
        } else {
            r_color = 1; g_color = 1; b_color = 1;
        }
    } else {
        // 单选模式或多选模式：支持相位显示
        const phaseOn = window.ElectronCloud.state.usePhaseColoring;
        let sign = 0;

        if (phaseOn) {
            let psi = 0;
            if (isMultiselectMode) {
                // 多选模式：只计算当前轨道的相位
                const R = Hydrogen.radialR(p.n, p.l, r);
                const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                psi = R * Y;
            } else {
                // 单选模式：计算所有轨道的叠加相位
                for (const params of paramsList) {
                    const R = Hydrogen.radialR(params.n, params.l, r);
                    const Y = Hydrogen.realYlm_value(params.angKey.l, params.angKey.m, params.angKey.t, theta, phi);
                    psi += R * Y;
                }
            }
            sign = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
        }

        // 使用全局 phaseColors 常量
        const pc = window.ElectronCloud.constants.phaseColors;
        if (sign > 0) {
            r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
        } else if (sign < 0) {
            r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
        } else {
            r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
        }
    }

    colors[i3] = r_color;
    colors[i3 + 1] = g_color;
    colors[i3 + 2] = b_color;

    // 同步更新 baseColors（用于星空闪烁模式）
    if (state.baseColors) {
        state.baseColors[i3] = r_color;
        state.baseColors[i3 + 1] = g_color;
        state.baseColors[i3 + 2] = b_color;
        state.baseColorsCount = state.pointCount + 1;
    }

    // 存储轨道索引
    if (state.pointOrbitalIndices) {
        state.pointOrbitalIndices[state.pointCount] = isIndependentMode ? orbitalIndex : 0;
    }

    // 记录轨道点映射（用于比照模式开关）
    if (isIndependentMode) {
        // 【关键修复】比照模式下构建完整键名，与 rolling update 和 data panel 保持一致
        const isCompareMode = window.ElectronCloud.ui.compareToggle && window.ElectronCloud.ui.compareToggle.checked;
        let orbitalKey;

        if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
            const slot = state.compareMode.activeSlots[orbitalIndex];
            // 必需格式: atom_orbital_slotIndex
            orbitalKey = slot ? `${slot.atom || 'H'}_${slot.orbital}_slot${slot.slotIndex}` : state.currentOrbitals[orbitalIndex];
        } else {
            orbitalKey = state.currentOrbitals[orbitalIndex];
        }

        if (orbitalKey) {
            if (!state.orbitalPointsMap[orbitalKey]) {
                state.orbitalPointsMap[orbitalKey] = [];
            }
            state.orbitalPointsMap[orbitalKey].push(state.pointCount);

            // 记录采样数据
            if (!state.orbitalSamplesMap[orbitalKey]) {
                state.orbitalSamplesMap[orbitalKey] = [];
            }
            state.orbitalSamplesMap[orbitalKey].push({
                r: r,
                theta: theta,
                probability: 0
            });
        }
    }

    state.pointCount++;
    state.radialSamples.push(r);
    state.angularSamples.push(theta);
    state.phiSamples.push(phi);

    window.ElectronCloud.Sampling.updateFarthestDistance(r);

    return true;
};

// ==================== 杂化模式采样 ====================

// 缓存杂化轨道的最大密度估计值
let hybridMaxDensityCache = null;
let hybridParamsHashCache = null;

/**
 * 计算参数列表的哈希值（用于缓存验证）
 */
function computeParamsHash(paramsList) {
    return paramsList.map(p => `${p.n}_${p.l}_${p.angKey.m}_${p.angKey.t}`).join('|');
}

/**
 * 杂化模式采样处理函数
 * 
 * 【物理原理】
 * 杂化轨道是多个原子轨道波函数的线性组合：Ψ_hybrid = Σ c_i × Ψ_i
 * 概率密度为 |Ψ_hybrid|² = |Σ c_i × R_i(r) × Y_i(θ,φ)|²
 * 
 * 【采样方法】
 * 优先使用高效的混合提议分布采样（hybridPreciseSample）
 * - 基于杂化轨道的径向CDF采样
 * - 然后使用角向权重接受-拒绝
 * 
 * 如果高效方法不可用，回退到基础拒绝采样
 * 
 * @param {Array} paramsList - 轨道参数列表
 * @param {Float32Array} positions - 位置数组
 * @param {Float32Array} colors - 颜色数组
 * @returns {boolean} - 是否成功添加了一个点
 */
window.ElectronCloud.Sampling.processHybridModePoint = function (paramsList, positions, colors) {
    const state = window.ElectronCloud.state;

    // 【关键修复】对轨道进行排序，确保符合杂化系数矩阵的预期顺序
    // 必须使用副本，不要修改原始 paramsList
    // 注意：sortOrbitalsForHybridization 返回的是新数组
    const sortedParamsList = Hydrogen.sortOrbitalsForHybridization ?
        Hydrogen.sortOrbitalsForHybridization(paramsList) : paramsList;

    // 使用排序后的列表进行后续计算
    paramsList = sortedParamsList;

    // 【新逻辑】使用 Round-Robin 轮流采样每条杂化轨道（均分点数）
    const numOrbitals = paramsList.length;

    // 初始化状态（如果需要）
    if (state.hybridRoundRobinIndex === undefined) {
        state.hybridRoundRobinIndex = 0;
    }
    if (!state.hybridOrbitalPointsMap) {
        state.hybridOrbitalPointsMap = {};
    }

    // 轮流选择杂化轨道
    const k = state.hybridRoundRobinIndex % numOrbitals;
    state.hybridRoundRobinIndex++;

    const coeffMatrix = Hydrogen.getHybridCoefficients(paramsList);
    const coeffs = coeffMatrix[k];

    // 将系数注入 paramsList
    const paramsWithCoeffs = paramsList.map((p, i) => ({
        ...p,
        coefficient: coeffs[i]
    }));

    // 采样
    // 获取当前原子类型
    const currentAtom = state.currentAtom || 'H';
    const result = Hydrogen.hybridPreciseSample(paramsWithCoeffs, state.samplingBoundary, currentAtom);

    if (result && result.accepted) {
        // 采样成功，记录这个点属于哪条杂化轨道
        result.hybridOrbitalIndex = k;
        return addHybridPointAll(result, paramsList, positions, colors);
    }

    // 如果高效采样失败（极少情况），回退到基础拒绝采样
    return processHybridPointFallback(paramsList, positions, colors);
};

/**
 * 处理单个杂化轨道的采样
 */
function processSingleHybridPoint(paramsList, hybridIndex, positions, colors) {
    const state = window.ElectronCloud.state;

    // 1. 获取杂化系数
    const numOrbitals = paramsList.length;
    const coeffMatrix = Hydrogen.getHybridCoefficients(paramsList);
    const coeffs = coeffMatrix[hybridIndex % numOrbitals];

    // 2. 将系数注入 paramsList
    const paramsWithCoeffs = paramsList.map((p, i) => ({
        ...p,
        coefficient: coeffs[i]
    }));

    // 3. 使用CDF采样
    // 这比基础拒绝采样快得多，且数值上更稳定
    const currentAtom = state.currentAtom || 'H';
    const result = Hydrogen.hybridPreciseSample(paramsWithCoeffs, state.samplingBoundary, currentAtom);

    if (result && result.accepted) {
        return addHybridPoint(result, paramsList, positions, colors);
    }

    return false;
}

/**
 * 估计单个杂化轨道的最大密度
 */
function estimateSingleHybridMaxDensity(paramsList, hybridIndex, rMax, numSamples = 2000) {
    let maxDensity = 0;

    for (let i = 0; i < numSamples; i++) {
        const r = Math.random() * rMax * Math.pow(Math.random(), 0.5);
        const cosTheta = 2 * Math.random() - 1;
        const theta = Math.acos(cosTheta);
        const phi = 2 * Math.PI * Math.random();

        const currentAtom = window.ElectronCloud.state.currentAtom || 'H';
        const density = Hydrogen.singleHybridDensity3D(paramsList, hybridIndex, r, theta, phi, 1, 1, currentAtom);
        if (density > maxDensity) {
            maxDensity = density;
        }
    }

    return maxDensity * 1.5;
}

/**
 * 添加杂化采样点到场景
 */
function addHybridPoint(result, paramsList, positions, colors) {
    const state = window.ElectronCloud.state;
    const i3 = state.pointCount * 3;

    positions[i3] = result.x;
    positions[i3 + 1] = result.y;
    positions[i3 + 2] = result.z;

    // 计算颜色（支持相位显示）
    let r_color, g_color, b_color;
    const phaseOn = state.usePhaseColoring;

    if (phaseOn) {
        // 使用已计算的波函数值确定相位
        const psi = result.psi;

        // 使用全局 phaseColors 常量
        const pc = window.ElectronCloud.constants.phaseColors;
        if (psi > 0) {
            r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
        } else if (psi < 0) {
            r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
        } else {
            r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
        }
    } else {
        // 默认白色
        r_color = 1; g_color = 1; b_color = 1;
    }

    colors[i3] = r_color;
    colors[i3 + 1] = g_color;
    colors[i3 + 2] = b_color;

    // 同步更新 baseColors（用于星空闪烁模式）
    if (state.baseColors) {
        state.baseColors[i3] = r_color;
        state.baseColors[i3 + 1] = g_color;
        state.baseColors[i3 + 2] = b_color;
        state.baseColorsCount = state.pointCount + 1;
    }

    // 存储轨道索引（杂化模式为 -1，表示混合轨道）
    if (state.pointOrbitalIndices) {
        state.pointOrbitalIndices[state.pointCount] = -1;
    }

    state.pointCount++;
    state.radialSamples.push(result.r);
    state.angularSamples.push(result.theta);
    state.phiSamples.push(result.phi);

    window.ElectronCloud.Sampling.updateFarthestDistance(result.r);

    return true;
}

/**
 * 回退方案：基础拒绝采样 - 【全部杂化轨道模式】
 * 
 * 计算所有 N 条杂化轨道的总概率密度：Σ|ψ_hybrid_i|²
 */
function processHybridPointFallback(paramsList, positions, colors) {
    const state = window.ElectronCloud.state;

    // 检查缓存是否有效
    const currentHash = computeParamsHash(paramsList);
    if (hybridParamsHashCache !== currentHash) {
        // 重新估计最大密度（使用所有杂化轨道的总密度）
        const rMax = Hydrogen.hybridRecommendRmax(paramsList);
        hybridMaxDensityCache = estimateAllHybridsMaxDensity(paramsList, rMax, 2000);
        hybridParamsHashCache = currentHash;
        console.log('杂化模式（全部）：估计最大密度 =', hybridMaxDensityCache.toExponential(3));
    }

    const samplingVolumeSize = state.samplingBoundary * 2;

    // 均匀随机采样
    const x = (Math.random() - 0.5) * samplingVolumeSize;
    const y = (Math.random() - 0.5) * samplingVolumeSize;
    const z = (Math.random() - 0.5) * samplingVolumeSize;
    const r = Math.sqrt(x * x + y * y + z * z);

    if (r < 1e-10) return false;

    const theta = Math.acos(z / r);
    const phi = Math.atan2(y, x);

    // 【关键修复】使用所有杂化轨道的总概率密度
    const currentAtom = state.currentAtom || 'H';
    const density = Hydrogen.allHybridOrbitalsDensity3D(paramsList, r, theta, phi, 1, 1, currentAtom);

    // 拒绝采样：接受概率 = density / maxDensity
    const acceptProb = density / hybridMaxDensityCache;

    // 动态调整最大密度估计（如果发现更大的密度）
    if (acceptProb > 1.0) {
        hybridMaxDensityCache = density * 1.5;
        console.log('杂化模式（全部）：更新最大密度估计 =', hybridMaxDensityCache.toExponential(3));
    }

    if (Math.random() > acceptProb) {
        return false; // 拒绝
    }

    // 确定此点属于哪条杂化轨道（用于着色）
    // 计算每条杂化轨道在此点的密度，选择最大的那条
    const numOrbitals = paramsList.length;
    let maxOrbitalDensity = 0;
    let dominantOrbitalIndex = 0;
    let dominantPsi = 0;

    for (let i = 0; i < numOrbitals; i++) {
        const orbitalDensity = Hydrogen.singleHybridDensity3D(paramsList, i, r, theta, phi, 1, 1, currentAtom);
        if (orbitalDensity > maxOrbitalDensity) {
            maxOrbitalDensity = orbitalDensity;
            dominantOrbitalIndex = i;
            dominantPsi = Hydrogen.singleHybridWavefunction(paramsList, i, r, theta, phi, 1, 1, currentAtom);
        }
    }

    // 添加点
    const result = {
        x: x + (Math.random() - 0.5) * 0.01,
        y: y + (Math.random() - 0.5) * 0.01,
        z: z + (Math.random() - 0.5) * 0.01,
        r, theta, phi, psi: dominantPsi,
        hybridOrbitalIndex: dominantOrbitalIndex, // 记录属于哪条杂化轨道
        accepted: true
    };

    return addHybridPointAll(result, paramsList, positions, colors);
}

/**
 * 估计所有杂化轨道总密度的最大值
 */
function estimateAllHybridsMaxDensity(paramsList, rMax, numSamples = 2000) {
    let maxDensity = 0;

    for (let i = 0; i < numSamples; i++) {
        const r = Math.random() * rMax * Math.pow(Math.random(), 0.5);
        const cosTheta = 2 * Math.random() - 1;
        const theta = Math.acos(cosTheta);
        const phi = 2 * Math.PI * Math.random();

        const currentAtom = window.ElectronCloud.state.currentAtom || 'H';
        const density = Hydrogen.allHybridOrbitalsDensity3D(paramsList, r, theta, phi, 1, 1, currentAtom);
        if (density > maxDensity) {
            maxDensity = density;
        }
    }

    return maxDensity * 1.5;
}

/**
 * 添加杂化采样点到场景 - 【全部模式】
 * 根据主导轨道索引着色，便于区分不同方向的杂化轨道
 */
function addHybridPointAll(result, paramsList, positions, colors) {
    const state = window.ElectronCloud.state;
    const i3 = state.pointCount * 3;

    positions[i3] = result.x;
    positions[i3 + 1] = result.y;
    positions[i3 + 2] = result.z;

    // 计算颜色
    let r_color, g_color, b_color;
    const phaseOn = state.usePhaseColoring;

    if (phaseOn) {
        // 相位模式：根据波函数符号着色
        const psi = result.psi;
        // 使用全局 phaseColors 常量
        const pc = window.ElectronCloud.constants.phaseColors;
        if (psi > 0) {
            r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
        } else if (psi < 0) {
            r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
        } else {
            r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
        }
    } else {
        // 默认白色
        r_color = 1; g_color = 1; b_color = 1;
    }

    colors[i3] = r_color;
    colors[i3 + 1] = g_color;
    colors[i3 + 2] = b_color;

    if (state.baseColors) {
        state.baseColors[i3] = r_color;
        state.baseColors[i3 + 1] = g_color;
        state.baseColors[i3 + 2] = b_color;
        state.baseColorsCount = state.pointCount + 1;
    }

    if (state.pointOrbitalIndices) {
        state.pointOrbitalIndices[state.pointCount] = result.hybridOrbitalIndex !== undefined ? result.hybridOrbitalIndex : -1;
    }

    // 【新增】记录点与杂化轨道的对应关系（用于可见性控制）
    if (result.hybridOrbitalIndex !== undefined) {
        const orbitalIndex = result.hybridOrbitalIndex;
        if (!state.hybridOrbitalPointsMap) {
            state.hybridOrbitalPointsMap = {};
        }
        if (!state.hybridOrbitalPointsMap[orbitalIndex]) {
            state.hybridOrbitalPointsMap[orbitalIndex] = [];
        }
        state.hybridOrbitalPointsMap[orbitalIndex].push(state.pointCount);
    }

    state.pointCount++;
    state.radialSamples.push(result.r);
    state.angularSamples.push(result.theta);
    state.phiSamples.push(result.phi);

    window.ElectronCloud.Sampling.updateFarthestDistance(result.r);

    return true;
}

/**
 * 重置杂化模式的缓存
 * 在切换轨道选择时调用
 */
window.ElectronCloud.Sampling.resetHybridCache = function () {
    hybridMaxDensityCache = null;
    hybridParamsHashCache = null;
    // 同时重置单个杂化轨道的缓存
    const state = window.ElectronCloud.state;
    if (state) {
        state._singleHybridMaxDensity = null;
        // 【新增】重置 Round-Robin 索引和轨道点映射
        state.hybridRoundRobinIndex = 0;
        state.hybridOrbitalPointsMap = {};
    }
    console.log('杂化模式：缓存已重置');
};

// 滚动生成更新逻辑
window.ElectronCloud.Sampling.performRollingUpdate = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const constants = window.ElectronCloud.constants;

    // 基本检查
    if (!state.points || state.pointCount < 100) return;

    // 每帧更新千分之一的点
    const pointsToUpdate = Math.max(1, Math.floor(state.pointCount / 1000)); // 0.1%

    // 获取必要的数据引用
    const positions = state.points.geometry.attributes.position.array;
    const colors = state.points.geometry.attributes.color.array;

    // 确定当前模式
    const isHybridMode = state.isHybridMode;
    const isMultiselectMode = ui.multiselectToggle && ui.multiselectToggle.checked;
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    const isIndependentMode = !isHybridMode && (isCompareMode || isMultiselectMode);

    // 获取轨道参数列表
    // 获取轨道参数列表
    // 【修复】移除全局 filter(Boolean)，保持索引一致性（特别是对比模式）
    let paramsList = state.currentOrbitals.map(k => Hydrogen.orbitalParamsFromKey(k));

    // 杂化模式需要纯净列表
    if (isHybridMode) {
        paramsList = paramsList.filter(Boolean);
    }

    if (!paramsList.length) return;

    const now = performance.now();

    // 获取当前原子类型
    const currentAtom = state.currentAtom || 'H';

    // 杂化模式预处理
    let hybridCoeffs = null;
    let hybridNumOrbitals = 0;
    if (isHybridMode) {
        if (Hydrogen.sortOrbitalsForHybridization) {
            paramsList = Hydrogen.sortOrbitalsForHybridization(paramsList);
        }
        hybridNumOrbitals = paramsList.length;
        const matrix = Hydrogen.getHybridCoefficients(paramsList);
        // 我们只在循环里取需要的系数
        hybridCoeffs = matrix;
    }

    // 比照模式：预先计算可见 slot（activeSlots 的索引），并用轮询确保每个轨道都会被更新
    const compareActiveSlots = (isCompareMode && state.compareMode && Array.isArray(state.compareMode.activeSlots)) ? state.compareMode.activeSlots : [];
    const compareSlotVisibility = (isCompareMode && state.compareMode) ? (state.compareMode.slotVisibility || {}) : {};
    const compareVisibleSlotIndices = [];
    if (isCompareMode && compareActiveSlots.length > 0) {
        for (let idx = 0; idx < compareActiveSlots.length; idx++) {
            const slot = compareActiveSlots[idx];
            if (!slot) continue;
            if (compareSlotVisibility[slot.slotIndex] !== false) compareVisibleSlotIndices.push(idx);
        }
    }
    if (state.rollingMode && state.rollingMode.compareRoundRobinIndex === undefined) {
        state.rollingMode.compareRoundRobinIndex = 0;
    }

    // 比照模式：严格按轮询选择“被替换点”（同一 slot 内替换）
    function pickCompareTargetIndex(pickedOrbitalIndex, orbitalKey) {
        if (!state.pointOrbitalIndices || state.pointCount <= 0) return -1;

        // 1) 先从该轨道的点集合里抽样（O(1) 期望）
        const list = state.orbitalPointsMap ? state.orbitalPointsMap[orbitalKey] : null;
        if (Array.isArray(list) && list.length > 0) {
            const tries = Math.min(30, list.length);
            for (let t = 0; t < tries; t++) {
                const cand = list[Math.floor(Math.random() * list.length)];
                if (cand === undefined || cand === null) continue;
                if (cand < 0 || cand >= state.pointCount) continue;
                if (state.pointOrbitalIndices[cand] === pickedOrbitalIndex) return cand;
            }
        }

        // 2) 兜底：用游标做有限步扫描，避免“随机扫不到”导致看起来只更新一个轨道
        if (!state.rollingMode) state.rollingMode = {};
        let cursor = Number.isInteger(state.rollingMode.compareScanCursor)
            ? state.rollingMode.compareScanCursor
            : Math.floor(Math.random() * state.pointCount);

        const maxSteps = Math.min(state.pointCount, 2000);
        for (let step = 0; step < maxSteps; step++) {
            const idx = (cursor + step) % state.pointCount;
            if (state.pointOrbitalIndices[idx] === pickedOrbitalIndex) {
                state.rollingMode.compareScanCursor = (idx + 1) % state.pointCount;
                return idx;
            }
        }

        state.rollingMode.compareScanCursor = (cursor + maxSteps) % state.pointCount;
        return -1;
    }

    for (let k = 0; k < pointsToUpdate; k++) {
        // 1) 选择要被替换的点：
        // - 比照模式：严格轮询每个可见 slot，并且只在该 slot 内替换
        // - 其他模式：全局随机替换
        let forcedCompare = null;
        let targetIndex = -1;

        if (isCompareMode && compareVisibleSlotIndices.length > 0 && isIndependentMode) {
            const rr = (state.rollingMode && Number.isInteger(state.rollingMode.compareRoundRobinIndex)) ? state.rollingMode.compareRoundRobinIndex : 0;
            const pickedOrbitalIndex = compareVisibleSlotIndices[rr % compareVisibleSlotIndices.length];
            if (state.rollingMode) state.rollingMode.compareRoundRobinIndex = rr + 1;

            const slot = compareActiveSlots[pickedOrbitalIndex];
            if (!slot) continue;
            const key = `${slot.atom || 'H'}_${slot.orbital}_slot${slot.slotIndex}`;
            forcedCompare = { orbitalIndex: pickedOrbitalIndex, orbitalKey: key, atomType: slot.atom || 'H' };

            targetIndex = pickCompareTargetIndex(pickedOrbitalIndex, key);
            if (targetIndex < 0) continue;
        } else {
            targetIndex = Math.floor(Math.random() * state.pointCount);
        }

        // 【关键修复】在覆盖 pointOrbitalIndices 之前缓存旧索引与对应的样本 key
        // 用于正确执行 orbitalSamplesMap 的 "One In, One Out"，避免滚动生成时图表看似不更新。
        const prevOrbitalIdx = state.pointOrbitalIndices ? state.pointOrbitalIndices[targetIndex] : -1;
        let prevOrbitalKey = null;
        if (prevOrbitalIdx !== undefined && prevOrbitalIdx >= 0) {
            if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                const prevSlot = state.compareMode.activeSlots[prevOrbitalIdx];
                if (prevSlot) {
                    prevOrbitalKey = `${prevSlot.atom || 'H'}_${prevSlot.orbital}_slot${prevSlot.slotIndex}`;
                }
            } else if (state.currentOrbitals && prevOrbitalIdx < state.currentOrbitals.length) {
                prevOrbitalKey = state.currentOrbitals[prevOrbitalIdx];
            }
        }

        // 【关键修复】跳过属于隐藏轨道的点，防止隐藏轨道的点被"蚕食"
        if (state.pointOrbitalIndices) {
            const checkOrbitalIdx = state.pointOrbitalIndices[targetIndex];

            if (isHybridMode && state.hybridOrbitalVisibility) {
                // 杂化模式：检查该点所属的轨道是否被隐藏
                if (checkOrbitalIdx !== undefined &&
                    state.hybridOrbitalVisibility[checkOrbitalIdx] === false) {
                    continue; // 跳过，不替换隐藏轨道的点
                }
            } else if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                // 【比照模式修复】使用slotVisibility检查
                // checkOrbitalIdx 是 Worker/paramsList 中的索引，对应 activeSlots 的索引
                const slotConfig = state.compareMode.activeSlots[checkOrbitalIdx];
                if (slotConfig && state.compareMode.slotVisibility) {
                    if (state.compareMode.slotVisibility[slotConfig.slotIndex] === false) {
                        continue; // 跳过，不替换隐藏轨道的点
                    }
                }
            } else if (isIndependentMode && state.orbitalVisibility && state.currentOrbitals) {
                // 独立模式（多选）：检查该点所属的轨道是否被隐藏
                const checkKey = state.currentOrbitals[checkOrbitalIdx];
                if (checkKey && state.orbitalVisibility[checkKey] === false) {
                    continue; // 跳过，不替换隐藏轨道的点
                }
            }
        }

        // 【关键修复】更新映射表（在覆盖前移除旧映射）
        // 比照模式 + forcedCompare（同轨道内替换）：不应把点从 orbitalPointsMap 移来移去，否则会导致某些轨道“只出不进”/统计偏置
        const isCompareSameOrbitalReplace = !!(isCompareMode && forcedCompare && prevOrbitalKey && forcedCompare.orbitalKey && prevOrbitalKey === forcedCompare.orbitalKey);
        if (!isCompareSameOrbitalReplace && state.pointOrbitalIndices) {
            const oldOrbitalIdx = prevOrbitalIdx;

            if (isHybridMode && state.hybridOrbitalPointsMap) {
                // 杂化模式：从旧的杂化轨道map中移除
                // oldOrbitalIdx 直接是杂化轨道索引 (0, 1, ...)
                if (oldOrbitalIdx !== undefined && oldOrbitalIdx >= 0) {
                    const arr = state.hybridOrbitalPointsMap[oldOrbitalIdx];
                    if (arr) {
                        const idxInMap = arr.indexOf(targetIndex);
                        if (idxInMap !== -1) arr.splice(idxInMap, 1);
                    }
                }
            } else if (isIndependentMode && state.orbitalPointsMap) {
                // 普通独立模式（比照/多选）：从旧的轨道map中移除
                if (oldOrbitalIdx !== undefined && oldOrbitalIdx >= 0 && prevOrbitalKey) {
                    const oldMap = state.orbitalPointsMap[prevOrbitalKey];
                    if (oldMap) {
                        const idxInMap = oldMap.indexOf(targetIndex);
                        if (idxInMap !== -1) {
                            try {
                                oldMap.splice(idxInMap, 1);
                            } catch (e) {
                                // 忽略潜在错误
                            }
                        }
                    }
                }
            }
        }

        // 2. 准备采样
        let x, y, z, r, theta, phi;
        let orbitalIndex = 0;
        let orbitalKey = null;
        let currentCoeffs = null; // 杂化系数

        if (isHybridMode) {
            // 【杂化采样】
            // 【关键修复】只为可见的轨道生成新点
            // 先收集所有可见的轨道索引
            const visibleHybridOrbitals = [];
            for (let h = 0; h < hybridNumOrbitals; h++) {
                if (!state.hybridOrbitalVisibility || state.hybridOrbitalVisibility[h] !== false) {
                    visibleHybridOrbitals.push(h);
                }
            }

            // 如果没有可见的轨道，跳过此次采样
            if (visibleHybridOrbitals.length === 0) continue;

            // 从可见轨道中随机选择一个
            orbitalIndex = visibleHybridOrbitals[Math.floor(Math.random() * visibleHybridOrbitals.length)];
            currentCoeffs = hybridCoeffs[orbitalIndex];

            // 注入系数
            const paramsWithCoeffs = paramsList.map((p, i) => ({
                ...p,
                coefficient: currentCoeffs[i]
            }));

            const result = Hydrogen.hybridPreciseSample(paramsWithCoeffs, state.samplingBoundary, currentAtom);
            if (!result || !result.accepted) continue;
            ({ x, y, z, r, theta, phi } = result);
            orbitalKey = state.currentOrbital;
        } else {
            // 【普通采样】
            if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                // 【比照模式】只为可见的 slot 生成新点
                if (forcedCompare) {
                    orbitalIndex = forcedCompare.orbitalIndex;
                } else {
                    const visibleSlotIndices = [];
                    state.compareMode.activeSlots.forEach((slot, idx) => {
                        const slotVis = state.compareMode.slotVisibility || {};
                        if (slotVis[slot.slotIndex] !== false) {
                            visibleSlotIndices.push(idx); // 保存 activeSlots 的索引
                        }
                    });

                    if (visibleSlotIndices.length === 0) continue;
                    orbitalIndex = visibleSlotIndices[Math.floor(Math.random() * visibleSlotIndices.length)];
                }
            } else if (isIndependentMode) {
                // 【多选模式】独立模式下也只为可见轨道生成新点
                // 收集可见轨道
                const visibleOrbitals = [];
                for (let o = 0; o < paramsList.length; o++) {
                    const key = state.currentOrbitals[o];
                    if (!state.orbitalVisibility || state.orbitalVisibility[key] !== false) {
                        visibleOrbitals.push(o);
                    }
                }

                // 如果没有可见轨道，跳过
                if (visibleOrbitals.length === 0) continue;

                orbitalIndex = visibleOrbitals[Math.floor(Math.random() * visibleOrbitals.length)];
            } else {
                orbitalIndex = 0;
            }

            const p = paramsList[orbitalIndex];
            if (!p) continue; // 【修复】防止 undefined 参数导致的 crash

            // 【关键修复】比照模式下从 activeSlots 获取 orbitalKey，确保与 handleWorkerResult 完全一致
            if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                if (forcedCompare && forcedCompare.orbitalKey) {
                    orbitalKey = forcedCompare.orbitalKey;
                } else {
                    const slot = state.compareMode.activeSlots[orbitalIndex];
                    // 必须构建完整键名：atom_orbital_slotIndex
                    orbitalKey = slot ? `${slot.atom || 'H'}_${slot.orbital}_slot${slot.slotIndex}` : state.currentOrbitals[orbitalIndex];
                }
            } else {
                orbitalKey = state.currentOrbitals[orbitalIndex];
            }

            // 【关键修复】比照模式下使用正确的原子类型
            let atomType = currentAtom;
            if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
                if (forcedCompare && forcedCompare.atomType) {
                    atomType = forcedCompare.atomType;
                } else {
                    const slot = state.compareMode.activeSlots[orbitalIndex];
                    if (slot) atomType = slot.atom || 'H';
                }
            }

            // 比照模式：同一轨道内替换时，尽量多尝试几次以保证每个轨道都能持续更新
            let result = null;
            const maxTries = (isCompareMode && forcedCompare) ? 12 : 1;
            for (let t = 0; t < maxTries; t++) {
                result = Hydrogen.importanceSample(p.n, p.l, p.angKey, state.samplingBoundary, atomType);
                if (result && result.accepted) break;
            }
            if (!result || !result.accepted) continue;
            ({ x, y, z, r, theta, phi } = result);
        }

        if (state.debugRollingCounter === undefined) state.debugRollingCounter = 0;
        state.debugRollingCounter++;

        // 3. 更新位置
        const i3 = targetIndex * 3;
        positions[i3] = x;
        positions[i3 + 1] = y;
        positions[i3 + 2] = z;

        // 【关键修复】同步更新 radialSamples 和 angularSamples
        // 确保直方图和能量分析图表能反映滚动更新后的数据！
        // targetIndex 是被替换的点的索引，直接覆盖
        if (state.radialSamples) state.radialSamples[targetIndex] = r;
        if (state.angularSamples) state.angularSamples[targetIndex] = theta;
        if (state.phiSamples) state.phiSamples[targetIndex] = phi;

        // 【关键修复】同步更新 originalPositions（用于位置隐藏方式的恢复）
        if (state.originalPositions) {
            state.originalPositions[i3] = x;
            state.originalPositions[i3 + 1] = y;
            state.originalPositions[i3 + 2] = z;
        }

        // 【关键修复】更新新映射
        if (state.pointOrbitalIndices) {
            // 杂化模式下 orbitalIndex 是 hybridIndex
            // 独立模式下 orbitalIndex 是 paramsList index
            const newIndexToStore = isHybridMode ? orbitalIndex : (isIndependentMode ? orbitalIndex : 0);
            state.pointOrbitalIndices[targetIndex] = newIndexToStore;
        }

        if (isIndependentMode && state.orbitalPointsMap && orbitalKey) {
            // 比照模式同轨道内替换：不需要修改 orbitalPointsMap（点不应在轨道之间漂移）
            if (!isCompareSameOrbitalReplace) {
                if (!state.orbitalPointsMap[orbitalKey]) state.orbitalPointsMap[orbitalKey] = [];
                state.orbitalPointsMap[orbitalKey].push(targetIndex);
            }

            // 【v11.0 修复】滚动生成模式下需要更新 orbitalSamplesMap 以驱动图表
            if (state.orbitalSamplesMap) {
                const indexMap = state.orbitalSampleIndexByPoint;
                const entry = indexMap ? indexMap[targetIndex] : null;

                // 比照模式同轨道内替换：就地更新对应样本（最稳定，不引入删加偏置）
                if (isCompareSameOrbitalReplace) {
                    const arr = state.orbitalSamplesMap[orbitalKey];
                    if (arr) {
                        // 若映射缺失/失效，先按 pointIndex 扫描一次来修复映射
                        if ((!entry || entry.key !== orbitalKey || !Number.isInteger(entry.idx)) && state.orbitalSampleIndexByPoint) {
                            for (let i = 0; i < arr.length; i++) {
                                const s = arr[i];
                                if (s && s.pointIndex === targetIndex) {
                                    state.orbitalSampleIndexByPoint[targetIndex] = { key: orbitalKey, idx: i };
                                    break;
                                }
                            }
                        }

                        const fixedEntry = state.orbitalSampleIndexByPoint ? state.orbitalSampleIndexByPoint[targetIndex] : null;
                        if (fixedEntry && fixedEntry.key === orbitalKey && Number.isInteger(fixedEntry.idx) && fixedEntry.idx >= 0 && fixedEntry.idx < arr.length) {
                            const s = arr[fixedEntry.idx];
                            if (s) {
                                s.r = r;
                                s.theta = theta;
                                s.phi = phi;
                                s.pointIndex = targetIndex;
                            } else {
                                arr[fixedEntry.idx] = { r, theta, phi, probability: 0, pointIndex: targetIndex };
                            }
                        } else {
                            // 极端回退：仍然保持集合规模稳定，追加并写回映射
                            arr.push({ r, theta, phi, probability: 0, pointIndex: targetIndex });
                            if (!state.orbitalSampleIndexByPoint) state.orbitalSampleIndexByPoint = {};
                            state.orbitalSampleIndexByPoint[targetIndex] = { key: orbitalKey, idx: arr.length - 1 };
                        }
                    }
                } else {
                    // 1. 移除旧样本（如果存在且 key 不同）
                    if (prevOrbitalKey && state.orbitalSamplesMap[prevOrbitalKey]) {
                        const oldSamples = state.orbitalSamplesMap[prevOrbitalKey];
                        const removeEntry = indexMap ? indexMap[targetIndex] : null;

                        if (removeEntry && removeEntry.key === prevOrbitalKey && Number.isInteger(removeEntry.idx) && removeEntry.idx >= 0 && removeEntry.idx < oldSamples.length) {
                            const removeIdx = removeEntry.idx;
                            const lastIdx = oldSamples.length - 1;
                            if (removeIdx !== lastIdx) {
                                const swapped = oldSamples[lastIdx];
                                oldSamples[removeIdx] = swapped;
                                if (swapped && swapped.pointIndex !== undefined && indexMap) {
                                    indexMap[swapped.pointIndex] = { key: prevOrbitalKey, idx: removeIdx };
                                }
                            }
                            oldSamples.pop();
                            if (indexMap) delete indexMap[targetIndex];
                        }
                    }

                    // 2. 添加新样本
                    if (orbitalKey) {
                        if (!state.orbitalSamplesMap[orbitalKey]) state.orbitalSamplesMap[orbitalKey] = [];
                        const arr = state.orbitalSamplesMap[orbitalKey];
                        arr.push({ r, theta, phi, probability: 0, pointIndex: targetIndex });
                        if (!state.orbitalSampleIndexByPoint) state.orbitalSampleIndexByPoint = {};
                        state.orbitalSampleIndexByPoint[targetIndex] = { key: orbitalKey, idx: arr.length - 1 };
                    }
                }
            }
        }

        if (isHybridMode && state.hybridOrbitalPointsMap) {
            if (!state.hybridOrbitalPointsMap[orbitalIndex]) state.hybridOrbitalPointsMap[orbitalIndex] = [];
            state.hybridOrbitalPointsMap[orbitalIndex].push(targetIndex);
            // 杂化模式通常重绘整个分布，或者暂不支持细粒度滚动更新图表
        }

        // 4. 计算颜色
        let r_color, g_color, b_color;

        if (isIndependentMode && !isMultiselectMode) {
            // 比照模式
            // 【关键修复】颜色绑定到 UI slotIndex，避免空 slot 造成颜色前移
            let colorIndex = orbitalIndex;
            if (isCompareMode && state.compareMode && Array.isArray(state.compareMode.activeSlots)) {
                const slot = state.compareMode.activeSlots[orbitalIndex];
                if (slot && Number.isInteger(slot.slotIndex)) colorIndex = slot.slotIndex;
            }

            if (constants.compareColors && constants.compareColors.length > 0) {
                const safe = ((colorIndex % constants.compareColors.length) + constants.compareColors.length) % constants.compareColors.length;
                const c = constants.compareColors[safe].value;
                r_color = c[0]; g_color = c[1]; b_color = c[2];
            } else {
                r_color = 1; g_color = 1; b_color = 1;
            }
        } else {
            // 相位颜色
            const phaseOn = state.usePhaseColoring;
            let sign = 0;

            if (phaseOn) {
                if (isHybridMode) {
                    // 杂化相位计算
                    let psi = 0;
                    for (let i = 0; i < paramsList.length; i++) {
                        const p = paramsList[i];
                        const coeff = currentCoeffs[i];
                        const R = Hydrogen.radialR(p.n, p.l, r, 1, 1, currentAtom);
                        const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        psi += coeff * R * Y;
                    }
                    sign = psi > 0 ? 1 : -1;
                } else if (isMultiselectMode) {
                    const p = paramsList[orbitalIndex];
                    const R = Hydrogen.radialR(p.n, p.l, r, 1, 1, currentAtom);
                    const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                    const psi = R * Y;
                    sign = psi > 0 ? 1 : -1;
                } else {
                    let psi = 0;
                    for (const p of paramsList) {
                        const R = Hydrogen.radialR(p.n, p.l, r, 1, 1, currentAtom);
                        const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                        psi += R * Y;
                    }
                    sign = psi > 0 ? 1 : -1;
                }
            }

            // 使用全局 phaseColors 常量，与所有采样函数保持一致
            const pc = constants.phaseColors;
            if (sign > 0) {
                r_color = pc.positive.r; g_color = pc.positive.g; b_color = pc.positive.b;
            } else if (sign < 0) {
                r_color = pc.negative.r; g_color = pc.negative.g; b_color = pc.negative.b;
            } else {
                r_color = pc.neutral.r; g_color = pc.neutral.g; b_color = pc.neutral.b;
            }
        }

        // 【关键修复】可见性检查
        // 如果轨道被隐藏，则强制颜色为不可见（黑/透明），且不触发高亮
        let isVisible = true;
        if (isHybridMode && state.hybridOrbitalVisibility) {
            // 杂化模式：检查杂化轨道索引对应的可见性
            if (state.hybridOrbitalVisibility[orbitalIndex] === false) {
                isVisible = false;
            }
        } else if (state.orbitalVisibility && orbitalKey) {
            // 普通模式：检查轨道key对应的可见性
            if (state.orbitalVisibility[orbitalKey] === false) {
                isVisible = false;
            }
        }

        if (!isVisible) {
            // 设为黑色 (配合 AdditiveBlending)
            r_color = 0; g_color = 0; b_color = 0;
        }

        // 同步 baseColors
        if (state.baseColors) {
            state.baseColors[i3] = r_color;
            state.baseColors[i3 + 1] = g_color;
            state.baseColors[i3 + 2] = b_color;
        }

        // 写入当前颜色（下一帧高亮前会短暂显示这个颜色，或者被高亮直接覆盖）
        colors[i3] = r_color;
        colors[i3 + 1] = g_color;
        colors[i3 + 2] = b_color;

        // 6. 注册高亮
        // 【关键修复】只有可见的点才高亮
        if (isVisible && state.rollingMode.highlightedPoints) {
            state.rollingMode.highlightedPoints.push({
                index: targetIndex,
                startTime: now
            });
        }

        // 7. 更新统计数据（确保图表能反映实时变化）
        if (state.radialSamples && targetIndex < state.radialSamples.length) {
            state.radialSamples[targetIndex] = r;
        }
        if (state.angularSamples && targetIndex < state.angularSamples.length) {
            state.angularSamples[targetIndex] = theta;
        }
        if (state.phiSamples && targetIndex < state.phiSamples.length) {
            state.phiSamples[targetIndex] = phi;
        }
    }

    // 【新增】图表动态更新 (仅在比照模式或多选模式且开启滚动生成时)
    if ((isCompareMode || isMultiselectMode) && state.rollingMode && state.rollingMode.enabled) {
        // 每200ms更新一次图表，与其他模式保持一致
        if (!state.lastChartUpdate || now - state.lastChartUpdate > 200) {
            // 【v11.0】orbitalSamplesMap 不再在滚动模式下累积，无需截断
            state.lastChartUpdate = now;
            if (window.ElectronCloud.Orbital.drawProbabilityChart) {
                // false 表示非 final 更新，可能触发轻量级渲染
                window.ElectronCloud.Orbital.drawProbabilityChart(false);
            }
        }
    }

    state.points.geometry.attributes.position.needsUpdate = true;
    state.points.geometry.attributes.color.needsUpdate = true;
};
