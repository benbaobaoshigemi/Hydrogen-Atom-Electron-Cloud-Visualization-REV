// 点采样逻辑模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Sampling = {};

// ==================== Web Worker 管理 ====================

// Worker 池配置
const WORKER_POOL_SIZE = navigator.hardwareConcurrency || 4; // 使用 CPU 核心数
const workerPool = [];
let workerReady = 0;
let pendingTasks = [];
let taskIdCounter = 0;
let useWorkers = true; // 是否使用 Worker（可降级）
let workersInitialized = false; // 防止重复初始化
let activeTasks = 0; // 当前正在执行的任务数
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
                    activeTasks--;

                    // 处理队列中的下一个任务
                    window.ElectronCloud.Sampling.processNextTask();
                }
            };

            worker.onerror = function (err) {
                console.error('Worker 错误:', err);
                worker.busy = false;
                activeTasks--;
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

    const positions = state.points.geometry.attributes.position.array;
    const colorsAttr = state.points.geometry.getAttribute('color');
    const colors = colorsAttr ? colorsAttr.array : null;

    // 确保 pointOrbitalIndices 数组存在
    if (!state.pointOrbitalIndices || state.pointOrbitalIndices.length < state.MAX_POINTS) {
        state.pointOrbitalIndices = new Int16Array(state.MAX_POINTS);
    }

    // 将 Worker 采样的点添加到主线程
    for (const point of result.points) {
        if (state.pointCount >= state.MAX_POINTS) break;

        const i3 = state.pointCount * 3;
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
        state.pointOrbitalIndices[state.pointCount] = point.orbitalIndex >= 0 ? point.orbitalIndex : 0;

        // 记录轨道点映射（用于比照模式开关）
        if (point.orbitalIndex >= 0) {
            const orbitalKey = state.currentOrbitals[point.orbitalIndex];
            if (orbitalKey) {
                if (!state.orbitalPointsMap[orbitalKey]) {
                    state.orbitalPointsMap[orbitalKey] = [];
                }
                state.orbitalPointsMap[orbitalKey].push(state.pointCount);
            }
        }

        state.pointCount++;
    }

    // 添加采样数据
    for (const sample of result.samples) {
        state.radialSamples.push(sample.r);
        state.angularSamples.push(sample.theta);

        if (sample.orbitalKey) {
            if (!state.orbitalSamplesMap[sample.orbitalKey]) {
                state.orbitalSamplesMap[sample.orbitalKey] = [];
            }
            state.orbitalSamplesMap[sample.orbitalKey].push({
                r: sample.r,
                theta: sample.theta,
                probability: 0 // Worker 中已计算过
            });
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

// 处理任务队列中的下一个任务
window.ElectronCloud.Sampling.processNextTask = function () {
    if (pendingTasks.length === 0) return;

    // 找到空闲的 Worker
    const availableWorker = workerPool.find(w => !w.busy);
    if (!availableWorker) return;

    const task = pendingTasks.shift();
    availableWorker.busy = true;
    activeTasks++;
    availableWorker.postMessage(task);
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
    const phaseOn = document.getElementById('phase-toggle')?.checked || false;

    // 准备比照颜色（只传递值，不传递对象）
    const compareColors = constants.compareColors ?
        constants.compareColors.map(c => c.value) : [];

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
                isHybridMode: isHybridMode, // 【新增】传递杂化模式标志
                phaseOn: phaseOn,
                compareColors: compareColors
            }
        };

        worker.busy = true;
        activeTasks++;
        worker.postMessage(task);
    }
};

// 终止所有 Worker
window.ElectronCloud.Sampling.terminateWorkers = function () {
    workerPool.forEach(w => w.terminate());
    workerPool.length = 0;
    workerReady = 0;
    pendingTasks = [];
    workersInitialized = false;
    activeTasks = 0;
};

// 重置采样会话（使旧的 Worker 结果无效）
window.ElectronCloud.Sampling.resetSamplingSession = function () {
    samplingSessionId++;
    pendingTasks = []; // 清空待处理队列
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

// 检查并设置采样模式
window.ElectronCloud.Sampling.setImportanceSamplingMode = function (enabled) {
    useImportanceSampling = enabled;
    console.log(`采样模式: ${enabled ? '重要性采样' : '传统拒绝采样'}`);
};

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
    const density = Hydrogen.density3D_real(p.angKey, p.n, p.l, r, theta, phi);

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
            const phaseOn = document.getElementById('phase-toggle')?.checked;
            let sign = 0;
            if (phaseOn) {
                // 只计算当前选中轨道的波函数，不叠加其他轨道
                const R = Hydrogen.radialR(p.n, p.l, r);
                const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                const psi = R * Y;
                sign = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
            }

            if (sign > 0) { // 红色（正相位）
                r_color = 1; g_color = 0.2; b_color = 0.2;
            } else if (sign < 0) { // 蓝色（负相位）
                r_color = 0.2; g_color = 0.2; b_color = 1;
            } else { // 中性白
                r_color = 1; g_color = 1; b_color = 1;
            }
        } else {
            // 比照模式：使用固定颜色区分不同轨道
            if (orbitalIndex < constants.compareColors.length) {
                const color = constants.compareColors[orbitalIndex].value;
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
    for (let i = 0; i < paramsList.length; i++) {
        const p = paramsList[i];
        const density = Hydrogen.density3D_real(p.angKey, p.n, p.l, r, theta, phi);
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
        const phaseOn = document.getElementById('phase-toggle')?.checked;
        let sign = 0;
        if (phaseOn) {
            let psi = 0;
            for (const p of paramsList) {
                const R = Hydrogen.radialR(p.n, p.l, r);
                const Y = Hydrogen.realYlm_value(p.angKey.l, p.angKey.m, p.angKey.t, theta, phi);
                psi += R * Y;
            }
            sign = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
        }

        let r_color, g_color, b_color;
        if (sign > 0) { // 红色
            r_color = 1; g_color = 0.2; b_color = 0.2;
        } else if (sign < 0) { // 蓝色
            r_color = 0.2; g_color = 0.2; b_color = 1;
        } else { // 中性白
            r_color = 1; g_color = 1; b_color = 1;
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

    // 使用重要性采样生成点（针对当前轨道优化）
    const result = Hydrogen.importanceSample(p.n, p.l, p.angKey, state.samplingBoundary);

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
        if (orbitalIndex < constants.compareColors.length) {
            const color = constants.compareColors[orbitalIndex].value;
            r_color = color[0];
            g_color = color[1];
            b_color = color[2];
        } else {
            r_color = 1; g_color = 1; b_color = 1;
        }
    } else {
        // 单选模式或多选模式：支持相位显示
        const phaseOn = document.getElementById('phase-toggle')?.checked;
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

        if (sign > 0) {
            r_color = 1; g_color = 0.2; b_color = 0.2;
        } else if (sign < 0) {
            r_color = 0.2; g_color = 0.2; b_color = 1;
        } else {
            r_color = 1; g_color = 1; b_color = 1;
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
        const orbitalKey = state.currentOrbitals[orbitalIndex];
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
 * - 基于杂化轨道的精确径向CDF采样
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

    const coeffMatrix = Hydrogen.getHybridCoefficients(numOrbitals);
    const coeffs = coeffMatrix[k];

    // 将系数注入 paramsList
    const paramsWithCoeffs = paramsList.map((p, i) => ({
        ...p,
        coefficient: coeffs[i]
    }));

    // 采样
    const result = Hydrogen.hybridPreciseSample(paramsWithCoeffs, state.samplingBoundary);

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
    const coeffMatrix = Hydrogen.getHybridCoefficients(numOrbitals);
    const coeffs = coeffMatrix[hybridIndex % numOrbitals];

    // 2. 将系数注入 paramsList
    const paramsWithCoeffs = paramsList.map((p, i) => ({
        ...p,
        coefficient: coeffs[i]
    }));

    // 3. 使用精确采样
    // 这比基础拒绝采样快得多，且物理上更准确
    const result = Hydrogen.hybridPreciseSample(paramsWithCoeffs, state.samplingBoundary);

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

        const density = Hydrogen.singleHybridDensity3D(paramsList, hybridIndex, r, theta, phi);
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
    const phaseOn = document.getElementById('phase-toggle')?.checked;

    if (phaseOn) {
        // 使用已计算的波函数值确定相位
        const psi = result.psi;

        if (psi > 0) {
            r_color = 1; g_color = 0.2; b_color = 0.2; // 正相位：红色
        } else if (psi < 0) {
            r_color = 0.2; g_color = 0.2; b_color = 1; // 负相位：蓝色
        } else {
            r_color = 1; g_color = 1; b_color = 1; // 节面：白色
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
    const density = Hydrogen.allHybridOrbitalsDensity3D(paramsList, r, theta, phi);

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
        const orbitalDensity = Hydrogen.singleHybridDensity3D(paramsList, i, r, theta, phi);
        if (orbitalDensity > maxOrbitalDensity) {
            maxOrbitalDensity = orbitalDensity;
            dominantOrbitalIndex = i;
            dominantPsi = Hydrogen.singleHybridWavefunction(paramsList, i, r, theta, phi);
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

        const density = Hydrogen.allHybridOrbitalsDensity3D(paramsList, r, theta, phi);
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
    const phaseOn = document.getElementById('phase-toggle')?.checked;

    if (phaseOn) {
        // 相位模式：根据波函数符号着色
        const psi = result.psi;
        if (psi > 0) {
            r_color = 1; g_color = 0.2; b_color = 0.2;
        } else if (psi < 0) {
            r_color = 0.2; g_color = 0.2; b_color = 1;
        } else {
            r_color = 1; g_color = 1; b_color = 1;
        }
    } else {
        // 默认白色（可以考虑未来根据 hybridOrbitalIndex 分配不同颜色）
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
