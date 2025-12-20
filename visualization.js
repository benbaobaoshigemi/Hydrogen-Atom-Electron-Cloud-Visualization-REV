// 3D 可视化效果模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Visualization = {};

// =============== 等值面计算 Worker 管理 ===============

let isosurfaceWorker = null;
let isosurfaceTaskId = 0;
const isosurfacePendingCallbacks = new Map();

// 初始化等值面 Worker
function initIsosurfaceWorker() {
    if (isosurfaceWorker) return;

    try {
        isosurfaceWorker = new Worker('isosurface-worker.js');

        isosurfaceWorker.onmessage = function (e) {
            const { type, taskId, progress, result, error } = e.data;

            if (type === 'ready') {
                console.log('等值面计算 Worker 就绪');
            } else if (type === 'progress') {
                // 进度回调
                const callbacks = isosurfacePendingCallbacks.get(taskId);
                if (callbacks && callbacks.onProgress) {
                    callbacks.onProgress(progress);
                }
            } else if (type === 'isosurface-result') {
                // 结果回调
                const callbacks = isosurfacePendingCallbacks.get(taskId);
                if (callbacks && callbacks.onComplete) {
                    callbacks.onComplete(result);
                }
                isosurfacePendingCallbacks.delete(taskId);
            } else if (type === 'error') {
                console.error('等值面计算 Worker 错误:', error);
                const callbacks = isosurfacePendingCallbacks.get(taskId);
                if (callbacks && callbacks.onError) {
                    callbacks.onError(error);
                }
                isosurfacePendingCallbacks.delete(taskId);
            }
        };

        isosurfaceWorker.onerror = function (err) {
            console.error('等值面计算 Worker 崩溃:', err);
        };

        console.log('等值面计算 Worker 初始化完成');
    } catch (err) {
        console.warn('等值面计算 Worker 初始化失败，将使用主线程计算:', err);
        isosurfaceWorker = null;
    }
}

// 异步计算等值面
window.ElectronCloud.Visualization.computeIsosurfaceAsync = function (params, onProgress, onComplete, onError) {
    const taskId = ++isosurfaceTaskId;

    if (!isosurfaceWorker) {
        initIsosurfaceWorker();
    }

    if (!isosurfaceWorker) {
        // Worker 不可用，同步计算
        if (onError) onError('Worker not available');
        return;
    }

    isosurfacePendingCallbacks.set(taskId, { onProgress, onComplete, onError });

    isosurfaceWorker.postMessage({
        type: 'compute-isosurface',
        taskId,
        data: params
    });

    return taskId;
};

// 检查 Worker 是否可用
window.ElectronCloud.Visualization.isIsosurfaceWorkerAvailable = function () {
    if (!isosurfaceWorker) initIsosurfaceWorker();
    return isosurfaceWorker !== null;
};

// 创建角向分布3D可视化
window.ElectronCloud.Visualization.createAngularOverlay = function () {
    const state = window.ElectronCloud.state;

    const overlayGroup = new THREE.Group();

    // 检查是否有选中的轨道
    if (!state.currentOrbitals || state.currentOrbitals.length === 0) return overlayGroup;

    // 【优化】角向显示大小设置为轨道轮廓（95%分位半径）的91%（原130%的70%）
    // 先计算95%分位半径
    let contourRadius = Math.max(15, state.farthestDistance * 0.6);
    if (window.ElectronCloud.Visualization.calculate95PercentileRadius) {
        const r95 = window.ElectronCloud.Visualization.calculate95PercentileRadius();
        if (r95 > 0) {
            contourRadius = r95;
        }
    }
    const baseRadius = contourRadius * 0.91; // 91% of contour (70% of original 130%)

    window.ElectronCloud.Visualization.createOrbitalMesh(overlayGroup, null, baseRadius);

    return overlayGroup;
};

/**
 * 创建轨道网格表面
 * 
 * 【关键】网格极点对准轨道的最高次旋转对称轴
 */
window.ElectronCloud.Visualization.createOrbitalMesh = function (group, params, baseRadius) {
    const state = window.ElectronCloud.state;
    const isHybridMode = state.isHybridMode === true;
    const isSingleHybridMode = isHybridMode && state.hybridRenderMode === 'single';
    const hybridIndex = state.hybridSingleIndex || 0;

    // 获取所有轨道参数
    let orbitalParams = state.currentOrbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
    // 【关键修复】对轨道进行排序，确保与 getHybridCoefficients 假设的顺序一致
    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }
    const numOrbitals = orbitalParams.length;

    // 获取杂化系数矩阵（如果是杂化模式）
    let coeffMatrix = null;
    if (isHybridMode && numOrbitals > 1 && Hydrogen.getHybridCoefficients) {
        coeffMatrix = Hydrogen.getHybridCoefficients(numOrbitals, orbitalParams);
    }

    // 【关键】确定网格极点的目标方向（最高次对称轴）
    let symmetryAxis = { x: 0, y: 0, z: 1 };

    if (isSingleHybridMode) {
        // 单个杂化轨道模式：使用解析算法计算当前 Lobe 的指向
        if (Hydrogen.getHybridLobeAxis) {
            symmetryAxis = Hydrogen.getHybridLobeAxis(orbitalParams, hybridIndex);
        } else if (Hydrogen.findHybridPrincipalAxis) {
            symmetryAxis = Hydrogen.findHybridPrincipalAxis(orbitalParams, hybridIndex);
        }
    } else if (isHybridMode && numOrbitals > 1) {
        // 全部杂化轨道模式：使用群论/成分分析计算对称轴
        if (Hydrogen.getSetSymmetryAxis) {
            symmetryAxis = Hydrogen.getSetSymmetryAxis(orbitalParams);
        } else if (Hydrogen.findHybridPrincipalAxis) {
            symmetryAxis = Hydrogen.findHybridPrincipalAxis(orbitalParams, 0);
        }
    } else if (!isHybridMode && numOrbitals === 1 && Hydrogen.getOrbitalSymmetryAxis) {
        // 单个原子轨道模式
        symmetryAxis = Hydrogen.getOrbitalSymmetryAxis(orbitalParams[0].angKey);
    }

    // 【关键】THREE.js SphereGeometry 极点在 Y 轴，需要旋转到对称轴方向
    let quaternion = null;
    const threeJsPoleAxis = new THREE.Vector3(0, 1, 0);
    const targetAxis = new THREE.Vector3(symmetryAxis.x, symmetryAxis.y, symmetryAxis.z).normalize();

    const dotProduct = threeJsPoleAxis.dot(targetAxis);
    if (Math.abs(dotProduct - 1) > 1e-6) {
        quaternion = new THREE.Quaternion();
        quaternion.setFromUnitVectors(threeJsPoleAxis, targetAxis);
    }

    // 角向强度计算函数（在物理坐标系中计算）
    function calcAngularIntensity(theta, phi) {
        if (isSingleHybridMode && coeffMatrix) {
            const coeffs = coeffMatrix[hybridIndex % numOrbitals];
            let psiSum = 0;
            for (let i = 0; i < numOrbitals; i++) {
                const op = orbitalParams[i];
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psiSum += coeffs[i] * Y;
            }
            return psiSum * psiSum; // |Ψ|²
        } else if (isHybridMode && numOrbitals > 1 && coeffMatrix) {
            let totalIntensity = 0;
            for (let h = 0; h < numOrbitals; h++) {
                const coeffs = coeffMatrix[h];
                let psiSum = 0;
                for (let i = 0; i < numOrbitals; i++) {
                    const op = orbitalParams[i];
                    const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                    psiSum += coeffs[i] * Y;
                }
                totalIntensity += psiSum * psiSum;
            }
            return totalIntensity;
        } else {
            let intensity = 0;
            for (const op of orbitalParams) {
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                intensity += Y * Y;
            }
            return intensity;
        }
    }

    // 提高网格密度 - 使用 IcosahedronGeometry (Icosphere) 替代 SphereGeometry
    // detail=5 (约 20,480 面) - 足够精细且性能良好
    const geometry = new THREE.IcosahedronGeometry(baseRadius, 50);
    const vertices = geometry.attributes.position.array;
    const colors = new Float32Array(vertices.length);

    for (let i = 0; i < vertices.length; i += 3) {
        const x0 = vertices[i];
        const y0 = vertices[i + 1];
        const z0 = vertices[i + 2];

        const r0 = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
        if (r0 < 1e-10) continue;

        // 将顶点方向旋转到物理坐标系
        let physDir = new THREE.Vector3(x0 / r0, y0 / r0, z0 / r0);
        if (quaternion) {
            physDir.applyQuaternion(quaternion);
        }

        const theta = Math.acos(Math.max(-1, Math.min(1, physDir.z)));
        const phi = Math.atan2(physDir.y, physDir.x);

        const totalAngularIntensity = calcAngularIntensity(theta, phi);

        const radiusScale = 1 + Math.sqrt(totalAngularIntensity) * 1.5;
        const newR = r0 * radiusScale;

        vertices[i] = physDir.x * newR;
        vertices[i + 1] = physDir.y * newR;
        vertices[i + 2] = physDir.z * newR;

        const intensity = Math.min(totalAngularIntensity * 3, 1);
        colors[i] = intensity;
        colors[i + 1] = intensity;
        colors[i + 2] = intensity;
    }

    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.attributes.position.needsUpdate = true;

    const material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.7,
        side: THREE.DoubleSide,
        wireframe: true,
        wireframeLinewidth: 2,
        depthWrite: false,
        blending: THREE.AdditiveBlending
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.layers.set(0);  // layer 0 以参与 Bloom 后处理
    group.add(mesh);
};

// 更新角向分布叠加
window.ElectronCloud.Visualization.updateAngularOverlay = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!ui.angular3dToggle || !ui.angular3dToggle.checked) return;

    if (state.angularOverlay) {
        state.scene.remove(state.angularOverlay);
        state.angularOverlay.traverse((child) => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
    }

    state.angularOverlay = window.ElectronCloud.Visualization.createAngularOverlay();
    state.angularOverlay.visible = true;

    if (state.points) {
        state.angularOverlay.rotation.copy(state.points.rotation);
        state.angularOverlay.updateMatrix();
    }

    state.scene.add(state.angularOverlay);
};

// 基于实际采样数据更新角向分布
window.ElectronCloud.Visualization.updateAngularOverlayFromSamples = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!ui.angular3dToggle || !ui.angular3dToggle.checked || !state.angularSamples || state.angularSamples.length === 0) return;

    if (state.angularOverlay) {
        state.scene.remove(state.angularOverlay);
        state.angularOverlay.traverse((child) => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
    }

    state.angularOverlay = window.ElectronCloud.Visualization.createAngularOverlayFromSamples();
    state.angularOverlay.visible = true;

    if (state.points) {
        state.angularOverlay.rotation.copy(state.points.rotation);
        state.angularOverlay.updateMatrix();
    }

    state.scene.add(state.angularOverlay);
};

// 基于采样数据创建角向分布
window.ElectronCloud.Visualization.createAngularOverlayFromSamples = function () {
    const state = window.ElectronCloud.state;

    const overlayGroup = new THREE.Group();

    if (!state.angularSamples || state.angularSamples.length === 0) {
        return window.ElectronCloud.Visualization.createAngularOverlay();
    }

    const baseRadius = Math.max(15, state.farthestDistance * 0.6);
    const thetaBins = 64;
    const phiBins = 128;

    const densityMap = new Array(thetaBins).fill(0).map(() => new Array(phiBins).fill(0));

    for (let i = 0; i < state.pointCount && i < state.radialSamples.length; i++) {
        const theta = state.angularSamples[i];

        if (state.points && state.points.geometry) {
            const positions = state.points.geometry.attributes.position.array;
            const x = positions[i * 3];
            const y = positions[i * 3 + 1];
            // const z = positions[i * 3 + 2]; 
            const phi = Math.atan2(y, x);

            const thetaIndex = Math.floor((theta / Math.PI) * (thetaBins - 1));
            const phiIndex = Math.floor(((phi + Math.PI) / (2 * Math.PI)) * (phiBins - 1));

            if (thetaIndex >= 0 && thetaIndex < thetaBins && phiIndex >= 0 && phiIndex < phiBins) {
                densityMap[thetaIndex][phiIndex]++;
            }
        }
    }

    let maxDensity = 0;
    for (let i = 0; i < thetaBins; i++) {
        for (let j = 0; j < phiBins; j++) {
            maxDensity = Math.max(maxDensity, densityMap[i][j]);
        }
    }

    if (maxDensity === 0) {
        return window.ElectronCloud.Visualization.createAngularOverlay();
    }

    // 使用 IcosahedronGeometry 替代 SphereGeometry，确保所有模式下网格密度一致且无极点
    // detail=5 (约 20,480 面) - 足够精细且性能良好
    const geometry = new THREE.IcosahedronGeometry(baseRadius, 50);
    const vertices = geometry.attributes.position.array;
    const colors = new Float32Array(vertices.length);

    for (let i = 0; i < vertices.length; i += 3) {
        const x = vertices[i];
        const y = vertices[i + 1];
        const z = vertices[i + 2];

        const r = Math.sqrt(x * x + y * y + z * z);
        const theta = Math.acos(z / r);
        const phi = Math.atan2(y, x);

        const thetaIndex = Math.floor((theta / Math.PI) * (thetaBins - 1));
        const phiIndex = Math.floor(((phi + Math.PI) / (2 * Math.PI)) * (phiBins - 1));

        let density = 0;
        if (thetaIndex >= 0 && thetaIndex < thetaBins && phiIndex >= 0 && phiIndex < phiBins) {
            density = densityMap[thetaIndex][phiIndex] / maxDensity;
        }

        const radiusScale = 1 + density * 2;
        vertices[i] = x * radiusScale;
        vertices[i + 1] = y * radiusScale;
        vertices[i + 2] = z * radiusScale;

        const intensityColor = Math.min(density * 1.2, 1);
        colors[i] = intensityColor;
        colors[i + 1] = intensityColor;
        colors[i + 2] = intensityColor;
    }

    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.attributes.position.needsUpdate = true;

    const material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.75,
        side: THREE.DoubleSide,
        wireframe: true,
        wireframeLinewidth: 2
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.layers.set(0);  // layer 0 以参与 Bloom 后处理
    overlayGroup.add(mesh);

    return overlayGroup;
};

// 导出截图
window.ElectronCloud.Visualization.exportImage = function () {
    const state = window.ElectronCloud.state;

    if (!state || !state.renderer) {
        alert('导出失败: 渲染器未初始化');
        return;
    }

    const backgroundSelect = document.getElementById('export-background-select');
    const background = backgroundSelect ? backgroundSelect.value : 'black';
    const canvas = state.renderer.domElement;

    try {
        if (background === 'transparent') {
            const tempCanvas = document.createElement('canvas');
            tempCanvas.width = canvas.width;
            tempCanvas.height = canvas.height;
            const ctx = tempCanvas.getContext('2d');
            ctx.drawImage(canvas, 0, 0);
            const imageData = ctx.getImageData(0, 0, tempCanvas.width, tempCanvas.height);
            const data = imageData.data;
            for (let i = 0; i < data.length; i += 4) {
                const r = data[i];
                const g = data[i + 1];
                const b = data[i + 2];
                const luminance = 0.299 * r + 0.587 * g + 0.114 * b;
                const maxChannel = Math.max(r, g, b);
                let alpha = Math.max(luminance, maxChannel);
                alpha = Math.min(255, alpha * 1.2);
                data[i + 3] = Math.round(alpha);
            }
            ctx.putImageData(imageData, 0, 0);
            tempCanvas.toBlob(function (blob) {
                if (!blob) return;
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                link.download = `electron-cloud-${canvas.width}x${canvas.height}-transparent-${timestamp}.png`;
                link.href = url;
                document.body.appendChild(link);
                link.click();
                setTimeout(() => { document.body.removeChild(link); URL.revokeObjectURL(url); }, 100);
            }, 'image/png');
        } else {
            const tempCanvas = document.createElement('canvas');
            tempCanvas.width = canvas.width;
            tempCanvas.height = canvas.height;
            const ctx = tempCanvas.getContext('2d');
            ctx.fillStyle = '#000000';
            ctx.fillRect(0, 0, tempCanvas.width, tempCanvas.height);
            ctx.drawImage(canvas, 0, 0);
            tempCanvas.toBlob(function (blob) {
                if (!blob) return;
                const url = URL.createObjectURL(blob);
                const link = document.createElement('a');
                const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
                link.download = `electron-cloud-${canvas.width}x${canvas.height}-black-${timestamp}.png`;
                link.href = url;
                document.body.appendChild(link);
                link.click();
                setTimeout(() => { document.body.removeChild(link); URL.revokeObjectURL(url); }, 100);
            }, 'image/png');
        }
    } catch (err) {
        console.error('导出截图失败:', err);
    }
};

// ==================== 95% 等值面轮廓（点云高亮方案）====================
// 【新方案】不使用网格，而是高亮现有采样点中靠近等值面的点
// 优势：天然解决波瓣分离问题，每个点独立，不会粘连

/**
 * 启用等值面轮廓高亮
 * 高亮现有点云中位于等值面附近的点
 */
window.ElectronCloud.Visualization.enableContourHighlight = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!state.points || state.pointCount < 100) {
        return;
    }

    // 标记等值面高亮已启用
    state.contourHighlightEnabled = true;

    // 计算每个点的密度排名（如果wave模式没有预先计算）
    const Hydrogen = window.Hydrogen;
    const positions = state.points.geometry.attributes.position;
    const colors = state.points.geometry.attributes.color;
    const pointCount = state.pointCount;
    const maxPoints = state.MAX_POINTS || 100000;

    // 获取轨道参数
    const isHybridMode = state.isHybridMode === true;
    const orbitals = state.currentOrbitals || [state.currentOrbital || '1s'];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    if (orbitalParams.length === 0) {
        return;
    }

    // 初始化等值面相关数组
    if (!state.contourRanks || state.contourRanks.length < maxPoints) {
        state.contourRanks = new Float32Array(maxPoints);
        state.contourDensities = new Float32Array(maxPoints);
        state.contourPsiSigns = new Int8Array(maxPoints); // 存储波函数符号
    }

    // 保存原始颜色（用于恢复）
    if (!state.contourBaseColors || state.contourBaseColors.length < maxPoints * 3) {
        state.contourBaseColors = new Float32Array(maxPoints * 3);
    }
    const colorArray = colors.array;
    for (let i = 0; i < pointCount * 3; i++) {
        state.contourBaseColors[i] = colorArray[i];
    }

    // 获取杂化系数矩阵
    let coeffMatrix = null;
    if (isHybridMode && Hydrogen.getHybridCoefficients && orbitalParams.length > 1) {
        coeffMatrix = Hydrogen.getHybridCoefficients(orbitalParams.length, orbitalParams);
    }

    // 是否为比照模式
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    const usePerPointOrbital = (isCompareMode || isHybridMode) && state.pointOrbitalIndices;

    // 计算每个点的密度和波函数符号（相位）
    const atomType = state.currentAtom || 'H';
    for (let i = 0; i < pointCount; i++) {
        const i3 = i * 3;
        const x = positions.array[i3];
        const y = positions.array[i3 + 1];
        const z = positions.array[i3 + 2];
        const r = Math.sqrt(x * x + y * y + z * z);
        const theta = r > 0 ? Math.acos(Math.max(-1, Math.min(1, z / r))) : 0;
        const phi = Math.atan2(y, x);

        let density = 0;
        let psi = 0;

        if (isHybridMode && coeffMatrix && usePerPointOrbital) {
            // 杂化模式：使用该点所属的杂化轨道
            const hybridIdx = state.pointOrbitalIndices[i] || 0;
            const coeffs = coeffMatrix[hybridIdx] || coeffMatrix[0];

            for (let j = 0; j < orbitalParams.length; j++) {
                const op = orbitalParams[j];
                const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += coeffs[j] * R * Y;
            }
            density = psi * psi;
        } else if (isCompareMode && usePerPointOrbital) {
            // 比照模式：使用该点所属的轨道
            const orbitalIndex = state.pointOrbitalIndices[i] || 0;
            const op = orbitalParams[orbitalIndex] || orbitalParams[0];
            const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
            const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
            psi = R * Y;
            density = psi * psi;
        } else {
            // 单选或多选叠加模式
            for (const op of orbitalParams) {
                const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += R * Y;
            }
            density = psi * psi;
        }

        state.contourDensities[i] = density;
        state.contourPsiSigns[i] = psi > 0 ? 1 : (psi < 0 ? -1 : 0);
    }

    // 按密度排序，计算每个点的等值面排名
    // 对于比照模式和杂化模式，按轨道分组排序
    if (usePerPointOrbital) {
        // 按轨道分组
        const orbitalPointGroups = {};
        for (let i = 0; i < pointCount; i++) {
            const orbitalIndex = state.pointOrbitalIndices[i] || 0;
            if (!orbitalPointGroups[orbitalIndex]) {
                orbitalPointGroups[orbitalIndex] = [];
            }
            orbitalPointGroups[orbitalIndex].push(i);
        }

        // 每组内部排序
        for (const groupKey of Object.keys(orbitalPointGroups)) {
            const groupIndices = orbitalPointGroups[groupKey];
            groupIndices.sort((a, b) => state.contourDensities[b] - state.contourDensities[a]);

            const groupSize = groupIndices.length;
            for (let rank = 0; rank < groupSize; rank++) {
                const pointIndex = groupIndices[rank];
                state.contourRanks[pointIndex] = rank / (groupSize - 1 || 1);
            }
        }
    } else {
        // 所有点一起排序
        const indices = new Array(pointCount);
        for (let i = 0; i < pointCount; i++) indices[i] = i;
        indices.sort((a, b) => state.contourDensities[b] - state.contourDensities[a]);

        for (let rank = 0; rank < pointCount; rank++) {
            const pointIndex = indices[rank];
            state.contourRanks[pointIndex] = rank / (pointCount - 1 || 1);
        }
    }

    // 应用等值面高亮效果
    window.ElectronCloud.Visualization.updateContourHighlight();
};

/**
 * 更新等值面高亮效果
 * 高亮排名在95%附近的点，并程序化插入更多点增加密度
 */
window.ElectronCloud.Visualization.updateContourHighlight = function () {
    const state = window.ElectronCloud.state;

    if (!state.contourHighlightEnabled || !state.points || !state.contourRanks) {
        return;
    }

    const colors = state.points.geometry.attributes.color;
    const positions = state.points.geometry.attributes.position;
    const pointCount = state.pointCount;
    const colorArray = colors.array;
    const posArray = positions.array;

    // 等值面位置：95%分位（rank = 0.95）
    const contourPosition = 0.95;
    // 等值面宽度：控制高亮区域的粗细（缩小到1%）
    const contourWidth = 0.01;

    // 高亮颜色（亮绿色）
    const highlightR = 0.3;
    const highlightG = 1.0;
    const highlightB = 0.6;

    // 收集等值面附近的点用于后续插点
    const contourPoints = [];

    for (let i = 0; i < pointCount; i++) {
        const i3 = i * 3;
        const rank = state.contourRanks[i];

        // 计算距离等值面的距离
        const distToContour = Math.abs(rank - contourPosition);

        if (distToContour < contourWidth) {
            // 在等值面附近：高亮显示
            // 使用更高的亮度（3.0-5.0）
            const intensity = 1.0 - (distToContour / contourWidth);
            const brightness = 3.0 + intensity * 2.0;

            colorArray[i3] = Math.min(1.0, highlightR * brightness);
            colorArray[i3 + 1] = Math.min(1.0, highlightG * brightness);
            colorArray[i3 + 2] = Math.min(1.0, highlightB * brightness);

            // 记录等值面点位置
            contourPoints.push({
                x: posArray[i3],
                y: posArray[i3 + 1],
                z: posArray[i3 + 2],
                psiSign: state.contourPsiSigns[i]
            });
        } else {
            // 远离等值面：显示为非常暗（几乎透明）
            const dimFactor = 0.05;
            colorArray[i3] = state.contourBaseColors[i3] * dimFactor;
            colorArray[i3 + 1] = state.contourBaseColors[i3 + 1] * dimFactor;
            colorArray[i3 + 2] = state.contourBaseColors[i3 + 2] * dimFactor;
        }
    }

    colors.needsUpdate = true;

    // 程序化插入更多点（在等值面点之间插值）
    if (contourPoints.length > 100) {
        window.ElectronCloud.Visualization.createContourInterpolationPoints(contourPoints);
    }
};

/**
 * 程序化生成等值面点（使用波函数直接计算）
 * 在球面上均匀采样方向，用二分搜索找到等值面半径，生成均匀分布的点
 */
window.ElectronCloud.Visualization.createContourInterpolationPoints = function (contourPoints) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const Hydrogen = window.Hydrogen;

    // 移除旧的插值点
    if (state.contourInterpolationPoints) {
        state.scene.remove(state.contourInterpolationPoints);
        state.contourInterpolationPoints.geometry.dispose();
        state.contourInterpolationPoints.material.dispose();
        state.contourInterpolationPoints = null;
    }

    // 获取轨道参数
    const isHybridMode = state.isHybridMode === true;
    const orbitals = state.currentOrbitals || [state.currentOrbital || '1s'];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    if (orbitalParams.length === 0) return;

    // 获取杂化系数
    let coeffMatrix = null;
    if (isHybridMode && Hydrogen.getHybridCoefficients && orbitalParams.length > 1) {
        coeffMatrix = Hydrogen.getHybridCoefficients(orbitalParams.length, orbitalParams);
    }

    // 计算波函数值
    const atomType = state.currentAtom || 'H';
    const calcPsi = (r, theta, phi) => {
        if (isHybridMode && coeffMatrix) {
            // 杂化模式：取所有杂化轨道中绝对值最大的波函数
            let maxAbsPsi = 0;
            let dominantPsi = 0;
            for (let h = 0; h < orbitalParams.length; h++) {
                const coeffs = coeffMatrix[h];
                let psi = 0;
                for (let j = 0; j < orbitalParams.length; j++) {
                    const op = orbitalParams[j];
                    const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                    const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                    psi += coeffs[j] * R * Y;
                }
                if (Math.abs(psi) > maxAbsPsi) {
                    maxAbsPsi = Math.abs(psi);
                    dominantPsi = psi;
                }
            }
            return dominantPsi;
        } else {
            // 单选或多选模式
            let psi = 0;
            for (const op of orbitalParams) {
                const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += R * Y;
            }
            return psi;
        }
    };

    // 计算密度
    const calcDensity = (r, theta, phi) => {
        const psi = calcPsi(r, theta, phi);
        return psi * psi;
    };

    // 从现有采样点估算95%密度阈值
    const sortedDensities = [...state.contourDensities.slice(0, state.pointCount)].sort((a, b) => b - a);
    const index95 = Math.floor(sortedDensities.length * 0.95);
    const densityThreshold = sortedDensities[index95] || 0;

    // 估算平均半径
    const positions = state.points.geometry.attributes.position.array;
    let sumR = 0;
    for (let i = 0; i < state.pointCount; i++) {
        const i3 = i * 3;
        const x = positions[i3];
        const y = positions[i3 + 1];
        const z = positions[i3 + 2];
        sumR += Math.sqrt(x * x + y * y + z * z);
    }
    const avgRadius = sumR / state.pointCount;

    // 在球面上均匀采样方向，生成等值面点
    const interpolatedPositions = [];
    const interpolatedColors = [];

    // 高亮颜色
    const highlightR = 0.3;
    const highlightG = 1.0;
    const highlightB = 0.6;
    const brightness = 4.0;

    // 使用黄金螺旋在球面上均匀分布点
    const numPoints = 30000; // 生成30000个均匀分布的方向
    const goldenAngle = Math.PI * (3 - Math.sqrt(5)); // 黄金角

    for (let i = 0; i < numPoints; i++) {
        // 均匀分布的 y 坐标 (-1 到 1)
        const y = 1 - (i / (numPoints - 1)) * 2;
        const radiusFactor = Math.sqrt(1 - y * y);

        // 均匀分布的 phi 角度
        const phi = goldenAngle * i;

        // 转换为球坐标的 theta（极角，从 z 轴测量）
        const theta = Math.acos(Math.max(-1, Math.min(1, y)));

        // 在这个方向上用二分搜索找等值面半径
        let rMin = 0.1;
        let rMax = avgRadius * 2.5;

        // 检查 rMax 处的密度是否足够低
        const densityAtMax = calcDensity(rMax, theta, phi);
        if (densityAtMax > densityThreshold) {
            rMax = avgRadius * 3.5;
        }

        // 二分搜索找等值面
        let radius95 = avgRadius;
        for (let iter = 0; iter < 20; iter++) {
            const rMid = (rMin + rMax) / 2;
            const density = calcDensity(rMid, theta, phi);

            if (density > densityThreshold) {
                rMin = rMid;
            } else {
                rMax = rMid;
            }
            radius95 = (rMin + rMax) / 2;
        }

        // 检查该方向的波函数相位
        const psi = calcPsi(radius95, theta, phi);

        // 只在有有效波函数的区域生成点（排除节面附近）
        if (Math.abs(psi) > densityThreshold * 0.1 && radius95 > 0.5) {
            // 转换为笛卡尔坐标
            const x = radius95 * Math.sin(theta) * Math.cos(phi);
            const yPos = radius95 * Math.cos(theta);
            const z = radius95 * Math.sin(theta) * Math.sin(phi);

            interpolatedPositions.push(x, yPos, z);
            interpolatedColors.push(
                Math.min(1.0, highlightR * brightness),
                Math.min(1.0, highlightG * brightness),
                Math.min(1.0, highlightB * brightness)
            );
        }
    }

    if (interpolatedPositions.length === 0) return;

    // 创建点云几何体
    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.Float32BufferAttribute(interpolatedPositions, 3));
    geometry.setAttribute('color', new THREE.Float32BufferAttribute(interpolatedColors, 3));

    // 使用与原始点云相同的材质
    const sprite = window.ElectronCloud.Scene.generateCircleSprite();
    const material = new THREE.PointsMaterial({
        map: sprite,
        size: state.points.material.size * 1.2,
        vertexColors: true,
        transparent: true,
        opacity: 1.0,
        depthWrite: false,
        blending: THREE.AdditiveBlending
    });

    state.contourInterpolationPoints = new THREE.Points(geometry, material);
    state.contourInterpolationPoints.layers.set(0);
    state.scene.add(state.contourInterpolationPoints);
};

/**
 * 禁用等值面轮廓高亮
 * 恢复原始点云颜色，并移除插值点
 */
window.ElectronCloud.Visualization.disableContourHighlight = function () {
    const state = window.ElectronCloud.state;

    // 移除插值点
    if (state.contourInterpolationPoints) {
        state.scene.remove(state.contourInterpolationPoints);
        state.contourInterpolationPoints.geometry.dispose();
        state.contourInterpolationPoints.material.dispose();
        state.contourInterpolationPoints = null;
    }

    if (!state.points || !state.contourBaseColors) {
        state.contourHighlightEnabled = false;
        return;
    }

    state.contourHighlightEnabled = false;

    // 恢复原始颜色
    const colors = state.points.geometry.attributes.color;
    const pointCount = state.pointCount;
    const colorArray = colors.array;

    for (let i = 0; i < pointCount * 3; i++) {
        colorArray[i] = state.contourBaseColors[i];
    }

    colors.needsUpdate = true;
};

// ==================== 网格等值面 ====================

window.ElectronCloud.Visualization.calculate95PercentileRadius = function () {
    const state = window.ElectronCloud.state;
    if (!state.radialSamples || state.radialSamples.length < 100) {
        return state.farthestDistance * 0.8;
    }
    const sortedRadii = state.radialSamples.slice(0, state.pointCount).sort((a, b) => a - b);
    const index95 = Math.floor(sortedRadii.length * 0.95);
    return sortedRadii[index95];
};

/**
 * 创建等值面网格 - 基于 Marching Cubes 算法
 * 每个波瓣生成独立的封闭等值面，自动分离连通域
 */
// 更新点云相位颜色
window.ElectronCloud.Visualization.updatePointColors = function () {
    const state = window.ElectronCloud.state;
    if (!state.points) return;

    const positions = state.points.geometry.attributes.position.array;
    const colors = state.points.geometry.attributes.color.array;
    const pointCount = state.pointCount;

    // 预备参数
    const isHybridMode = state.isHybridMode === true;
    const orbitals = state.currentOrbitals || [];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    if (orbitalParams.length === 0) return;

    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    let coeffMatrix = null;
    if (isHybridMode && orbitalParams.length > 1 && Hydrogen.getHybridCoefficients) {
        coeffMatrix = Hydrogen.getHybridCoefficients(orbitalParams.length, orbitalParams);
    }

    // 使用全局 phaseColors 常量，与所有采样函数保持一致
    const pc = window.ElectronCloud.constants.phaseColors;
    const colorPos = { r: pc.positive.r, g: pc.positive.g, b: pc.positive.b }; // 正相位：蓝色
    const colorNeg = { r: pc.negative.r, g: pc.negative.g, b: pc.negative.b }; // 负相位：红色
    const colorWhite = { r: pc.neutral.r, g: pc.neutral.g, b: pc.neutral.b };
    const atomType = state.currentAtom || 'H';

    // 遍历所有点
    for (let i = 0; i < pointCount; i++) {
        const i3 = i * 3;
        const x = positions[i3];
        const y = positions[i3 + 1];
        const z = positions[i3 + 2];

        // 默认白色 (未开启相位显示时)
        if (!state.usePhaseColoring) {
            colors[i3] = colorWhite.r;
            colors[i3 + 1] = colorWhite.g;
            colors[i3 + 2] = colorWhite.b;
            continue;
        }

        // 计算 Psi (复用 createContourMesh 的逻辑)
        let psi = 0;
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r < 1e-10) {
            psi = 0;
        } else {
            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);

            if (isHybridMode && coeffMatrix) {
                let maxAbsPsi = 0;
                for (let h = 0; h < orbitalParams.length; h++) {
                    let hPsi = 0;
                    for (let k = 0; k < orbitalParams.length; k++) {
                        const op = orbitalParams[k];
                        hPsi += coeffMatrix[h][k] * Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType) *
                            Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                    }
                    if (Math.abs(hPsi) > maxAbsPsi) {
                        maxAbsPsi = Math.abs(hPsi);
                        psi = hPsi;
                    }
                }
            } else {
                // 单轨/多轨叠加
                // 如果是多轨模式，点云中的点实际上属于特定轨道
                // 但为了相位显示的连续性，我们可以计算叠加场的相位，或者区分轨道
                // 这里采用：计算叠加波函数 psi (简单叠加)
                for (const op of orbitalParams) {
                    psi += Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType) *
                        Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                }
            }
        }

        // 根据相位上色
        if (psi > 0) {
            colors[i3] = colorPos.r;
            colors[i3 + 1] = colorPos.g;
            colors[i3 + 2] = colorPos.b;
        } else if (psi < 0) {
            colors[i3] = colorNeg.r;
            colors[i3 + 1] = colorNeg.g;
            colors[i3 + 2] = colorNeg.b;
        } else {
            colors[i3] = colorWhite.r;
            colors[i3 + 1] = colorWhite.g;
            colors[i3 + 2] = colorWhite.b;
        }
    }

    state.points.geometry.attributes.color.needsUpdate = true;
};

window.ElectronCloud.Visualization.createContourMesh = function (group, baseRadius) {
    const state = window.ElectronCloud.state;
    const isHybridMode = state.isHybridMode === true;

    const orbitals = state.currentOrbitals || [];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);
    if (orbitalParams.length === 0) return;

    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    let coeffMatrix = null;
    if (isHybridMode && orbitalParams.length > 1 && Hydrogen.getHybridCoefficients) {
        coeffMatrix = Hydrogen.getHybridCoefficients(orbitalParams.length, orbitalParams);
    }

    // 波函数计算 (直角坐标) - 用于计算 isovalue
    const atomType = state.currentAtom || 'H';
    function calcPsiLocal(x, y, z) {
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r < 1e-10) return 0;
        const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
        const phi = Math.atan2(y, x);

        if (isHybridMode && coeffMatrix) {
            let maxAbsPsi = 0, dominantPsi = 0;
            for (let h = 0; h < orbitalParams.length; h++) {
                let psi = 0;
                for (let i = 0; i < orbitalParams.length; i++) {
                    const op = orbitalParams[i];
                    psi += coeffMatrix[h][i] * Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType) *
                        Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                }
                if (Math.abs(psi) > maxAbsPsi) { maxAbsPsi = Math.abs(psi); dominantPsi = psi; }
            }
            return dominantPsi;
        } else {
            let psi = 0;
            for (const op of orbitalParams) {
                psi += Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType) *
                    Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
            }
            return psi;
        }
    }

    // 计算等值面阈值
    const psiValues = [];
    const positions = state.points.geometry.attributes.position.array;
    for (let i = 0; i < state.pointCount; i++) {
        const i3 = i * 3;
        psiValues.push(Math.abs(calcPsiLocal(positions[i3], positions[i3 + 1], positions[i3 + 2])));
    }
    psiValues.sort((a, b) => b - a);
    const isovalue = psiValues[Math.floor(psiValues.length * 0.95)] || 0.0001;

    const bound = baseRadius * 1.3;
    const resolution = 200; // 【提升】从 160 提升到 200

    // 检查 Worker 是否可用
    const useWorker = window.ElectronCloud.Visualization.isIsosurfaceWorkerAvailable &&
        window.ElectronCloud.Visualization.isIsosurfaceWorkerAvailable();

    if (useWorker) {
        // 【异步】使用 Worker 计算
        const atomSlaterData = (window.SlaterBasis && window.SlaterBasis[atomType])
            ? window.SlaterBasis[atomType] : null;

        const workerParams = {
            orbitalParams: orbitalParams.map(op => ({
                n: op.n, l: op.l,
                angKey: { l: op.angKey.l, m: op.angKey.m, t: op.angKey.t }
            })),
            coeffs: isHybridMode && coeffMatrix ? coeffMatrix[0] : orbitalParams.map(() => 1),
            bound: bound,
            resolution: resolution,
            isovalue: isovalue,
            atomType: atomType,
            slaterBasis: atomSlaterData,
            color: { r: 0.2, g: 1.0, b: 0.5 }
        };

        window.ElectronCloud.Visualization.computeIsosurfaceAsync(
            workerParams,
            null,
            function (result) {
                fillGroupWithWorkerResult(group, result, state.usePhaseColoring);
                console.log('普通模式等值面计算完成 (Worker)');
            },
            function (error) {
                console.error('Worker 等值面计算失败:', error);
            }
        );
    } else {
        // 【同步】回退到主线程计算
        const result = window.MarchingCubes.run(calcPsiLocal, { min: [-bound, -bound, -bound], max: [bound, bound, bound] }, resolution, isovalue);
        addLobeMeshesSync(group, result, state.usePhaseColoring);
    }

    // Worker 结果填充
    function fillGroupWithWorkerResult(grp, result, usePhase) {
        const colorPos = usePhase ? { r: 0.0, g: 0.75, b: 1.0 } : { r: 0.2, g: 1.0, b: 0.5 };
        const colorNeg = usePhase ? { r: 1.0, g: 0.27, b: 0.0 } : { r: 0.2, g: 1.0, b: 0.5 };

        for (const vertices of result.positive) {
            if (vertices.length < 9) continue;
            addMeshFromVertices(grp, vertices, colorPos);
        }
        for (const vertices of result.negative) {
            if (vertices.length < 9) continue;
            addMeshFromVertices(grp, vertices, colorNeg);
        }
    }

    function addMeshFromVertices(grp, vertices, color) {
        const posArr = new Float32Array(vertices);
        const geometry = new THREE.BufferGeometry();
        geometry.setAttribute('position', new THREE.BufferAttribute(posArr, 3));
        geometry.computeVertexNormals();

        const colors = new Float32Array(vertices.length);
        for (let i = 0; i < vertices.length / 3; i++) {
            colors[i * 3] = color.r;
            colors[i * 3 + 1] = color.g;
            colors[i * 3 + 2] = color.b;
        }
        geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

        const mesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({
            vertexColors: true, transparent: true, opacity: 0.55,
            side: THREE.DoubleSide, wireframe: true, wireframeLinewidth: 1.5
        }));
        mesh.layers.set(1);
        grp.add(mesh);
    }

    // 同步回退
    function addLobeMeshesSync(grp, result, usePhase) {
        const colorPos = usePhase ? { r: 0.0, g: 0.75, b: 1.0 } : { r: 0.2, g: 1.0, b: 0.5 };
        const colorNeg = usePhase ? { r: 1.0, g: 0.27, b: 0.0 } : { r: 0.2, g: 1.0, b: 0.5 };

        function addLobes(triangles, color) {
            if (triangles.length < 9) return;
            const components = window.MarchingCubes.separate(triangles);
            for (const comp of components) {
                if (comp.length < 9) continue;
                const geom = window.MarchingCubes.toGeometry(comp);
                const cols = new Float32Array(comp.length * 3);
                for (let i = 0; i < comp.length; i++) {
                    cols[i * 3] = color.r;
                    cols[i * 3 + 1] = color.g;
                    cols[i * 3 + 2] = color.b;
                }
                geom.setAttribute('color', new THREE.BufferAttribute(cols, 3));
                const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({
                    vertexColors: true, transparent: true, opacity: 0.55,
                    side: THREE.DoubleSide, wireframe: true, wireframeLinewidth: 1.5
                }));
                mesh.layers.set(1);
                grp.add(mesh);
            }
        }
        addLobes(result.positive, colorPos);
        addLobes(result.negative, colorNeg);
    }
};

window.ElectronCloud.Visualization.createContourOverlay = function () {
    const state = window.ElectronCloud.state;
    const group = new THREE.Group();
    const ui = window.ElectronCloud.ui;
    const constants = window.ElectronCloud.constants;

    if (!state.points) return group;

    // 计算包围球半径
    const positions = state.points.geometry.attributes.position.array;
    let maxR = 0;
    for (let i = 0; i < state.pointCount; i++) {
        const x = positions[i * 3];
        const y = positions[i * 3 + 1];
        const z = positions[i * 3 + 2];
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r > maxR) maxR = r;
    }
    const baseRadius = maxR > 0 ? maxR : 10;

    // 检查是否为比照模式
    const isCompareMode = ui && ui.compareToggle && ui.compareToggle.checked;

    if (isCompareMode && state.compareMode && state.compareMode.activeSlots) {
        // 比照模式：为每个slot创建独立的等值面，使用对应颜色
        window.ElectronCloud.Visualization.createCompareContourMeshes(group, baseRadius);
    } else {
        // 其他模式：使用原有逻辑
        window.ElectronCloud.Visualization.createContourMesh(group, baseRadius);
    }

    return group;
};

/**
 * 比照模式专用：为每个slot创建独立的等值面轮廓
 * 每个slot使用对应的颜色（红/绿/蓝）
 */
window.ElectronCloud.Visualization.createCompareContourMeshes = function (group, baseRadius) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    const Hydrogen = window.Hydrogen;

    const activeSlots = state.compareMode.activeSlots || [];
    if (activeSlots.length === 0) return;

    const compareColors = constants.compareColors || [
        { name: 'red', value: [1, 0.2, 0.2] },
        { name: 'green', value: [0.2, 1, 0.2] },
        { name: 'blue', value: [0.2, 0.2, 1] }
    ];

    // 为每个slot创建等值面
    for (let slotIdx = 0; slotIdx < activeSlots.length; slotIdx++) {
        const slot = activeSlots[slotIdx];
        if (!slot || !slot.orbital) continue;

        // 检查该slot是否可见
        if (state.compareMode.slotVisibility &&
            state.compareMode.slotVisibility[slot.slotIndex] === false) {
            continue;
        }

        const orbitalKey = slot.orbital;
        const atomType = slot.atom || 'H';
        const orbitalParams = Hydrogen.orbitalParamsFromKey(orbitalKey);

        if (!orbitalParams) continue;

        // 获取该slot的颜色
        const colorValue = compareColors[slotIdx % compareColors.length].value;
        const slotColor = { r: colorValue[0], g: colorValue[1], b: colorValue[2] };

        // 为该slot的波函数计算等值面
        function calcPsi(x, y, z) {
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r < 1e-10) return 0;
            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);

            return Hydrogen.radialWavefunction(orbitalParams.n, orbitalParams.l, r, 1, 1, atomType) *
                Hydrogen.realYlm_value(orbitalParams.angKey.l, orbitalParams.angKey.m,
                    orbitalParams.angKey.t, theta, phi);
        }

        // 收集该slot的采样点来计算isovalue
        const psiValues = [];
        const positions = state.points.geometry.attributes.position.array;

        // 使用pointOrbitalIndices找到属于该slot的点
        for (let i = 0; i < state.pointCount; i++) {
            // 检查该点是否属于当前slot
            if (state.pointOrbitalIndices && state.pointOrbitalIndices[i] === slotIdx) {
                const i3 = i * 3;
                const psi = calcPsi(positions[i3], positions[i3 + 1], positions[i3 + 2]);
                psiValues.push(Math.abs(psi));
            }
        }

        // 如果没有足够的点，使用全局采样估算
        if (psiValues.length < 50) {
            // 使用球面均匀采样
            const sampleCount = 1000;
            const radius = baseRadius * 0.7;
            for (let s = 0; s < sampleCount; s++) {
                const theta = Math.acos(2 * Math.random() - 1);
                const phi = 2 * Math.PI * Math.random();
                const r = Math.random() * radius;
                const x = r * Math.sin(theta) * Math.cos(phi);
                const y = r * Math.cos(theta);
                const z = r * Math.sin(theta) * Math.sin(phi);
                const psi = calcPsi(x, y, z);
                psiValues.push(Math.abs(psi));
            }
        }

        psiValues.sort((a, b) => b - a);
        const isovalue = psiValues[Math.floor(psiValues.length * 0.95)] || 0.0001;

        // Marching Cubes
        const bound = baseRadius * 1.3;
        const resolution = 160; // 与单轨模式统一

        const result = window.MarchingCubes.run(
            calcPsi,
            { min: [-bound, -bound, -bound], max: [bound, bound, bound] },
            resolution,
            isovalue
        );

        // 创建网格
        function addLobeMeshes(triangles, color) {
            if (triangles.length < 9) return;
            const components = window.MarchingCubes.separate(triangles);

            for (const comp of components) {
                if (comp.length < 9) continue;
                const geom = window.MarchingCubes.toGeometry(comp);
                const colors = new Float32Array(comp.length * 3);
                for (let i = 0; i < comp.length; i++) {
                    colors[i * 3] = color.r;
                    colors[i * 3 + 1] = color.g;
                    colors[i * 3 + 2] = color.b;
                }
                geom.setAttribute('color', new THREE.BufferAttribute(colors, 3));

                const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({
                    vertexColors: true, transparent: true, opacity: 0.55,
                    side: THREE.DoubleSide, wireframe: true, wireframeLinewidth: 1.5
                }));
                // 【关键新增】标记该网格属于哪个slot，以便后续单独控制可见性
                mesh.userData.slotIndex = slotIdx;
                mesh.layers.set(1);
                group.add(mesh);
            }
        }

        // 使用slot对应颜色
        addLobeMeshes(result.positive, slotColor);
        addLobeMeshes(result.negative, slotColor);
    }
};

window.ElectronCloud.Visualization.updateContourOverlay = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    const toggle = document.getElementById('contour-3d-toggle');
    if (!toggle || !toggle.checked) {
        if (state.contourOverlay) {
            state.scene.remove(state.contourOverlay);
            state.contourOverlay.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
            state.contourOverlay = null;
        }
        return;
    }

    if (state.contourOverlay) {
        state.scene.remove(state.contourOverlay);
        state.contourOverlay.traverse((child) => {
            if (child.geometry) child.geometry.dispose();
            if (child.material) child.material.dispose();
        });
    }

    state.contourOverlay = window.ElectronCloud.Visualization.createContourOverlay();
    state.contourOverlay.visible = true;

    if (state.points) {
        state.contourOverlay.rotation.copy(state.points.rotation);
        state.contourOverlay.updateMatrix();
    }

    state.scene.add(state.contourOverlay);
};

/**
 * 创建杂化轨道的独立等值面轮廓（每个轨道一个）
 * 【性能优化】使用 Worker 异步计算 Marching Cubes，避免阻塞主线程
 * 先返回占位 Group 对象，Worker 完成后自动填充几何体
 */
window.ElectronCloud.Visualization.createHybridContourOverlays = function () {
    const state = window.ElectronCloud.state;
    const overlays = [];

    if (!state.isHybridMode || !state.hybridOrbitalPointsMap) {
        return overlays;
    }

    const orbitals = state.currentOrbitals || [];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    if (orbitalParams.length === 0) return overlays;

    // 获取原子类型
    const atomType = state.currentAtom || 'H';

    // 对轨道进行排序
    if (Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    const numOrbitals = orbitalParams.length;
    const coeffMatrix = Hydrogen.getHybridCoefficients ? Hydrogen.getHybridCoefficients(numOrbitals, orbitalParams) : null;

    if (!coeffMatrix) return overlays;

    // 计算轨道边界
    let estimatedRadius = 15;
    if (Hydrogen.estimateOrbitalRadius95) {
        for (const key of orbitals) {
            const r95 = Hydrogen.estimateOrbitalRadius95(atomType, key);
            if (r95 > estimatedRadius) estimatedRadius = r95;
        }
    }

    // 颜色配置
    const colors = [
        { r: 0.2, g: 1.0, b: 0.5 },
        { r: 1.0, g: 0.5, b: 0.2 },
        { r: 0.2, g: 0.5, b: 1.0 },
        { r: 1.0, g: 0.2, b: 0.8 },
        { r: 0.8, g: 1.0, b: 0.2 },
        { r: 0.2, g: 1.0, b: 1.0 },
    ];

    // 计算统一的 isovalue
    let unifiedIsovalue = 0.0001;
    const positions = state.points.geometry.attributes.position.array;
    const allPsiValues = [];

    for (let hybridIndex = 0; hybridIndex < numOrbitals; hybridIndex++) {
        const pointIndices = state.hybridOrbitalPointsMap[hybridIndex] || [];
        const coeffs = coeffMatrix[hybridIndex];

        for (const pointIdx of pointIndices) {
            const i3 = pointIdx * 3;
            const x = positions[i3];
            const y = positions[i3 + 1];
            const z = positions[i3 + 2];
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r < 1e-10) continue;
            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);

            let psi = 0;
            for (let j = 0; j < orbitalParams.length; j++) {
                const op = orbitalParams[j];
                const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += coeffs[j] * R * Y;
            }
            allPsiValues.push(Math.abs(psi));
        }
    }

    if (allPsiValues.length >= 100) {
        allPsiValues.sort((a, b) => b - a);
        unifiedIsovalue = allPsiValues[Math.floor(allPsiValues.length * 0.95)] || 0.0001;
    }

    // 检查 Worker 是否可用
    const useWorker = window.ElectronCloud.Visualization.isIsosurfaceWorkerAvailable &&
        window.ElectronCloud.Visualization.isIsosurfaceWorkerAvailable();

    // 为每个杂化轨道创建占位 Group 或同步计算
    for (let hybridIndex = 0; hybridIndex < numOrbitals; hybridIndex++) {
        if (state.hybridOrbitalVisibility && state.hybridOrbitalVisibility[hybridIndex] === false) {
            overlays.push(null);
            continue;
        }

        const pointIndices = state.hybridOrbitalPointsMap[hybridIndex] || [];
        if (pointIndices.length < 50) {
            overlays.push(null);
            continue;
        }

        // 计算 maxR
        let maxR = 0;
        for (const pointIdx of pointIndices) {
            const i3 = pointIdx * 3;
            const x = positions[i3], y = positions[i3 + 1], z = positions[i3 + 2];
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r > maxR) maxR = r;
        }
        const bound = Math.max(maxR, estimatedRadius) * 1.5;
        const color = colors[hybridIndex % colors.length];
        const coeffs = coeffMatrix[hybridIndex];

        // 创建占位 Group
        const overlayGroup = new THREE.Group();
        overlayGroup.userData.hybridIndex = hybridIndex;
        overlays.push(overlayGroup);

        if (useWorker) {
            // 【异步】使用 Worker 计算
            // 获取该原子的 SlaterBasis 数据（如果有）
            const atomSlaterData = (window.SlaterBasis && window.SlaterBasis[atomType])
                ? window.SlaterBasis[atomType] : null;

            const workerParams = {
                orbitalParams: orbitalParams.map(op => ({
                    n: op.n, l: op.l,
                    angKey: { l: op.angKey.l, m: op.angKey.m, t: op.angKey.t }
                })),
                coeffs: coeffs,
                bound: bound,
                resolution: 200,
                isovalue: unifiedIsovalue,
                atomType: atomType,
                slaterBasis: atomSlaterData, // 传递 SlaterBasis 数据
                color: color
            };

            // 创建闭包捕获当前的 overlayGroup 和配置
            (function (group, hybridIdx, cfg) {
                window.ElectronCloud.Visualization.computeIsosurfaceAsync(
                    cfg,
                    null, // onProgress
                    function (result) {
                        // Worker 完成，填充几何体
                        fillGroupWithResult(group, result, cfg.color);
                        console.log(`杂化轨道 ${hybridIdx} 等值面计算完成 (Worker)`);
                    },
                    function (error) {
                        console.error(`Worker 等值面计算失败 (轨道 ${hybridIdx}):`, error);
                    }
                );
            })(overlayGroup, hybridIndex, workerParams);

        } else {
            // 【同步】回退到主线程计算
            fillGroupSync(overlayGroup, orbitalParams, coeffs, bound, unifiedIsovalue, atomType, color);
        }
    }

    return overlays;

    // 将 Worker 结果填充到 Group
    function fillGroupWithResult(group, result, color) {
        const usePhase = state.usePhaseColoring;
        const colorPos = usePhase ? { r: 0.0, g: 0.75, b: 1.0 } : color;
        const colorNeg = usePhase ? { r: 1.0, g: 0.27, b: 0.0 } : color;

        addComponentsToGroup(group, result.positive, colorPos);
        addComponentsToGroup(group, result.negative, colorNeg);
    }

    function addComponentsToGroup(group, components, meshColor) {
        for (const vertices of components) {
            if (vertices.length < 9) continue;

            const positions = new Float32Array(vertices);
            const geometry = new THREE.BufferGeometry();
            geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
            geometry.computeVertexNormals();

            const colors = new Float32Array(vertices.length);
            for (let i = 0; i < vertices.length / 3; i++) {
                colors[i * 3] = meshColor.r;
                colors[i * 3 + 1] = meshColor.g;
                colors[i * 3 + 2] = meshColor.b;
            }
            geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));

            const mesh = new THREE.Mesh(geometry, new THREE.MeshBasicMaterial({
                vertexColors: true, transparent: true, opacity: 0.6,
                side: THREE.DoubleSide, wireframe: true, wireframeLinewidth: 1.0
            }));
            mesh.layers.set(1);
            group.add(mesh);
        }
    }

    // 同步计算回退
    function fillGroupSync(group, orbitalParams, coeffs, bound, isovalue, atomType, color) {
        function calcPsi(x, y, z) {
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r < 1e-10) return 0;
            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);

            let psi = 0;
            for (let j = 0; j < orbitalParams.length; j++) {
                const op = orbitalParams[j];
                const R = Hydrogen.radialWavefunction(op.n, op.l, r, 1, 1, atomType);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += coeffs[j] * R * Y;
            }
            return psi;
        }

        const result = window.MarchingCubes.run(calcPsi, { min: [-bound, -bound, -bound], max: [bound, bound, bound] }, 160, isovalue);

        function addLobeMeshes(triangles, meshColor) {
            if (triangles.length < 9) return;
            const components = window.MarchingCubes.separate(triangles);

            for (const comp of components) {
                if (comp.length < 9) continue;
                const geom = window.MarchingCubes.toGeometry(comp);
                const colors = new Float32Array(comp.length * 3);
                for (let i = 0; i < comp.length; i++) {
                    colors[i * 3] = meshColor.r;
                    colors[i * 3 + 1] = meshColor.g;
                    colors[i * 3 + 2] = meshColor.b;
                }
                geom.setAttribute('color', new THREE.BufferAttribute(colors, 3));

                const mesh = new THREE.Mesh(geom, new THREE.MeshBasicMaterial({
                    vertexColors: true, transparent: true, opacity: 0.6,
                    side: THREE.DoubleSide, wireframe: true, wireframeLinewidth: 1.0
                }));
                mesh.layers.set(1);
                group.add(mesh);
            }
        }

        const usePhase = state.usePhaseColoring;
        const colorPos = usePhase ? { r: 0.0, g: 0.75, b: 1.0 } : color;
        const colorNeg = usePhase ? { r: 1.0, g: 0.27, b: 0.0 } : color;
        addLobeMeshes(result.positive, colorPos);
        addLobeMeshes(result.negative, colorNeg);
    }
};

/**
 * 更新比照模式下特定slot的等值面可见性
 * @param {number} slotIndex - 插槽索引
 * @param {boolean} visible - 是否可见
 */
window.ElectronCloud.Visualization.updateCompareContourVisibility = function (slotIndex, visible) {
    const state = window.ElectronCloud.state;
    // 检查轮廓层是否存在
    if (!state.contourOverlay) return;

    // 遍历所有子网格
    state.contourOverlay.traverse((child) => {
        // 检查是否有 slotIndex 标记
        if (child.userData && child.userData.slotIndex !== undefined) {
            if (child.userData.slotIndex === slotIndex) {
                child.visible = visible;
            }
        }
    });

    // 请求重绘
    state.sceneNeedsUpdate = true;
};

/**
 * 创建高密度测地线球体网格 (Icosphere)
 * 手动细分二十面体，绕过 Three.js 版本限制
 * @param { number } radius 半径
 * @param { number } detail 细分等级(建议 5 或 6)
 */
window.ElectronCloud.Visualization.createGeodesicIcosphere = function (radius, detail) {
    const t = (1 + Math.sqrt(5)) / 2;
    const vertices = [
        -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t, 0,
        0, -1, t, 0, 1, t, 0, -1, -t, 0, 1, -t,
        t, 0, -1, t, 0, 1, -t, 0, -1, -t, 0, 1
    ];
    const indices = [
        0, 11, 5, 0, 5, 1, 0, 1, 7, 0, 7, 10, 0, 10, 11,
        1, 5, 9, 5, 11, 4, 11, 10, 2, 10, 7, 6, 7, 1, 8,
        3, 9, 4, 3, 4, 2, 3, 2, 6, 3, 6, 8, 3, 8, 9,
        4, 9, 5, 2, 4, 11, 6, 2, 10, 8, 6, 7, 9, 8, 1
    ];

    let faces = [];
    for (let i = 0; i < indices.length; i += 3) {
        faces.push([
            new THREE.Vector3(vertices[indices[i] * 3], vertices[indices[i] * 3 + 1], vertices[indices[i] * 3 + 2]).normalize().multiplyScalar(radius),
            new THREE.Vector3(vertices[indices[i + 1] * 3], vertices[indices[i + 1] * 3 + 1], vertices[indices[i + 1] * 3 + 2]).normalize().multiplyScalar(radius),
            new THREE.Vector3(vertices[indices[i + 2] * 3], vertices[indices[i + 2] * 3 + 1], vertices[indices[i + 2] * 3 + 2]).normalize().multiplyScalar(radius)
        ]);
    }

    for (let d = 0; d < detail; d++) {
        const newFaces = [];
        for (const tri of faces) {
            const v1 = tri[0];
            const v2 = tri[1];
            const v3 = tri[2];
            const a = v1.clone().add(v2).multiplyScalar(0.5).normalize().multiplyScalar(radius);
            const b = v2.clone().add(v3).multiplyScalar(0.5).normalize().multiplyScalar(radius);
            const c = v3.clone().add(v1).multiplyScalar(0.5).normalize().multiplyScalar(radius);

            newFaces.push([v1, a, c]);
            newFaces.push([v2, b, a]);
            newFaces.push([v3, c, b]);
            newFaces.push([a, b, c]);
        }
        faces = newFaces;
    }

    const positions = new Float32Array(faces.length * 9);
    for (let i = 0; i < faces.length; i++) {
        positions[i * 9] = faces[i][0].x;
        positions[i * 9 + 1] = faces[i][0].y;
        positions[i * 9 + 2] = faces[i][0].z;
        positions[i * 9 + 3] = faces[i][1].x;
        positions[i * 9 + 4] = faces[i][1].y;
        positions[i * 9 + 5] = faces[i][1].z;
        positions[i * 9 + 6] = faces[i][2].x;
        positions[i * 9 + 7] = faces[i][2].y;
        positions[i * 9 + 8] = faces[i][2].z;
    }

    const geometry = new THREE.BufferGeometry();
    geometry.setAttribute('position', new THREE.BufferAttribute(positions, 3));
    return geometry;
};

