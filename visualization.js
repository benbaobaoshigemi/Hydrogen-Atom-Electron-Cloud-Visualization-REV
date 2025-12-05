// 3D 可视化效果模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.Visualization = {};

// 创建角向分布3D可视化
window.ElectronCloud.Visualization.createAngularOverlay = function () {
    const state = window.ElectronCloud.state;

    const overlayGroup = new THREE.Group();

    // 检查是否有选中的轨道
    if (!state.currentOrbitals || state.currentOrbitals.length === 0) return overlayGroup;

    // 【优化】角向显示大小设置为轨道轮廓（95%分位半径）的130%
    // 先计算95%分位半径
    let contourRadius = Math.max(15, state.farthestDistance * 0.6);
    if (window.ElectronCloud.Visualization.calculate95PercentileRadius) {
        const r95 = window.ElectronCloud.Visualization.calculate95PercentileRadius();
        if (r95 > 0) {
            contourRadius = r95;
        }
    }
    const baseRadius = contourRadius * 1.3; // 130% of contour

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
        coeffMatrix = Hydrogen.getHybridCoefficients(numOrbitals);
    }

    // 【关键】确定网格极点的目标方向（最高次对称轴）
    let symmetryAxis = { x: 0, y: 0, z: 1 };

    console.log('[DEBUG createOrbitalMesh] isHybridMode:', isHybridMode, 'isSingleHybridMode:', isSingleHybridMode, 'hybridIndex:', hybridIndex);

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

    // 提高网格密度
    const geometry = new THREE.SphereGeometry(baseRadius, 96, 48);
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
    mesh.layers.set(1);
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

    const geometry = new THREE.SphereGeometry(baseRadius, phiBins, thetaBins);
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
    mesh.layers.set(1);
    overlayGroup.add(mesh);

    return overlayGroup;
};

// 导出截图
window.ElectronCloud.Visualization.exportImage = function () {
    console.log('exportImage 函数被调用');
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

// ==================== 95% 等值面轮廓 ====================

window.ElectronCloud.Visualization.createContourOverlay = function () {
    const state = window.ElectronCloud.state;
    const overlayGroup = new THREE.Group();

    if (!state.points || state.pointCount < 100) {
        return overlayGroup;
    }

    const contourRadius = window.ElectronCloud.Visualization.calculate95PercentileRadius();

    if (contourRadius <= 0) {
        return overlayGroup;
    }

    window.ElectronCloud.Visualization.createContourMesh(overlayGroup, contourRadius);
    return overlayGroup;
};

window.ElectronCloud.Visualization.calculate95PercentileRadius = function () {
    const state = window.ElectronCloud.state;
    if (!state.radialSamples || state.radialSamples.length < 100) {
        return state.farthestDistance * 0.8;
    }
    const sortedRadii = state.radialSamples.slice(0, state.pointCount).sort((a, b) => a - b);
    const index95 = Math.floor(sortedRadii.length * 0.95);
    return sortedRadii[index95];
};

window.ElectronCloud.Visualization.calculate95PercentileRadiusMap = function (thetaBins = 32, phiBins = 64) {
    const state = window.ElectronCloud.state;
    const radiusBins = new Array(thetaBins).fill(null).map(() =>
        new Array(phiBins).fill(null).map(() => [])
    );

    if (!state.points || !state.points.geometry) {
        return null;
    }

    const positions = state.points.geometry.attributes.position.array;

    for (let i = 0; i < state.pointCount; i++) {
        const i3 = i * 3;
        const x = positions[i3];
        const y = positions[i3 + 1];
        const z = positions[i3 + 2];
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r < 1e-10) continue;
        const theta = Math.acos(z / r);
        const phi = Math.atan2(y, x);
        const thetaIndex = Math.min(thetaBins - 1, Math.floor((theta / Math.PI) * thetaBins));
        const phiIndex = Math.min(phiBins - 1, Math.floor(((phi + Math.PI) / (2 * Math.PI)) * phiBins));
        radiusBins[thetaIndex][phiIndex].push(r);
    }

    const radiusMap = new Array(thetaBins).fill(0).map(() => new Array(phiBins).fill(0));
    let globalMax = 0;

    for (let i = 0; i < thetaBins; i++) {
        for (let j = 0; j < phiBins; j++) {
            const samples = radiusBins[i][j];
            if (samples.length < 5) {
                radiusMap[i][j] = -1;
            } else {
                samples.sort((a, b) => a - b);
                const index95 = Math.floor(samples.length * 0.95);
                radiusMap[i][j] = samples[index95];
                globalMax = Math.max(globalMax, samples[index95]);
            }
        }
    }

    for (let i = 0; i < thetaBins; i++) {
        for (let j = 0; j < phiBins; j++) {
            if (radiusMap[i][j] < 0) {
                let sum = 0, count = 0;
                for (let di = -2; di <= 2; di++) {
                    for (let dj = -2; dj <= 2; dj++) {
                        const ni = (i + di + thetaBins) % thetaBins;
                        const nj = (j + dj + phiBins) % phiBins;
                        if (radiusMap[ni][nj] > 0) {
                            sum += radiusMap[ni][nj];
                            count++;
                        }
                    }
                }
                radiusMap[i][j] = count > 0 ? sum / count : globalMax * 0.8;
            }
        }
    }

    return { radiusMap, thetaBins, phiBins, globalMax };
};

/**
 * 创建等值面网格
 * 
 * 【关键】网格极点对准轨道的最高次旋转对称轴
 */
window.ElectronCloud.Visualization.createContourMesh = function (group, baseRadius) {
    const state = window.ElectronCloud.state;
    const isHybridMode = state.isHybridMode === true;
    const isSingleHybridMode = isHybridMode && state.hybridRenderMode === 'single';
    const hybridIndex = state.hybridSingleIndex || 0;

    // 获取所有轨道参数
    const orbitals = state.currentOrbitals || [];
    let orbitalParams = orbitals.map(key => Hydrogen.orbitalParamsFromKey(key)).filter(Boolean);

    if (orbitalParams.length === 0) {
        return;
    }

    // 【关键修复】对轨道进行排序
    if (isHybridMode && Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }
    const numOrbitals = orbitalParams.length;

    // 获取杂化系数矩阵
    let coeffMatrix = null;
    if (isHybridMode && numOrbitals > 1 && Hydrogen.getHybridCoefficients) {
        coeffMatrix = Hydrogen.getHybridCoefficients(numOrbitals);
    }

    // 【关键】确定网格极点的目标方向（最高次对称轴）
    let symmetryAxis = { x: 0, y: 0, z: 1 };

    if (isSingleHybridMode) {
        if (Hydrogen.getHybridLobeAxis) {
            symmetryAxis = Hydrogen.getHybridLobeAxis(orbitalParams, hybridIndex);
        } else if (Hydrogen.findHybridPrincipalAxis) {
            symmetryAxis = Hydrogen.findHybridPrincipalAxis(orbitalParams, hybridIndex);
        }
    } else if (isHybridMode && numOrbitals > 1) {
        if (Hydrogen.getSetSymmetryAxis) {
            symmetryAxis = Hydrogen.getSetSymmetryAxis(orbitalParams);
        } else if (Hydrogen.findHybridPrincipalAxis) {
            symmetryAxis = Hydrogen.findHybridPrincipalAxis(orbitalParams, 0);
        }
    } else if (!isHybridMode && numOrbitals === 1 && Hydrogen.getOrbitalSymmetryAxis) {
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

    // 密度计算函数
    function calcDensity(r, theta, phi) {
        if (isSingleHybridMode && Hydrogen.singleHybridDensity3D) {
            return Hydrogen.singleHybridDensity3D(orbitalParams, hybridIndex, r, theta, phi);
        } else if (isHybridMode && Hydrogen.allHybridOrbitalsDensity3D) {
            return Hydrogen.allHybridOrbitalsDensity3D(orbitalParams, r, theta, phi);
        } else {
            let density = 0;
            for (const op of orbitalParams) {
                density += Hydrogen.density3D_real(op.angKey, op.n, op.l, r, theta, phi);
            }
            return density / numOrbitals;
        }
    }

    // 计算95%分位的密度阈值
    const densities = [];
    const positions = state.points.geometry.attributes.position.array;

    for (let i = 0; i < state.pointCount; i++) {
        const i3 = i * 3;
        const x = positions[i3];
        const y = positions[i3 + 1];
        const z = positions[i3 + 2];
        const r = Math.sqrt(x * x + y * y + z * z);
        if (r < 1e-10) continue;

        const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
        const phi = Math.atan2(y, x);

        const density = calcDensity(r, theta, phi);
        densities.push(density);
    }

    densities.sort((a, b) => b - a);
    const index95 = Math.floor(densities.length * 0.95);
    const densityThreshold = densities[index95] || 0;

    const thetaBins = 48;
    const phiBins = 96;
    const geometry = new THREE.SphereGeometry(1, phiBins, thetaBins);
    const vertices = geometry.attributes.position.array;
    const colors = new Float32Array(vertices.length);

    for (let i = 0; i < vertices.length; i += 3) {
        const x0 = vertices[i];
        const y0 = vertices[i + 1];
        const z0 = vertices[i + 2];

        const unitR = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
        if (unitR < 1e-10) continue;

        let physDir = new THREE.Vector3(x0 / unitR, y0 / unitR, z0 / unitR);
        if (quaternion) {
            physDir.applyQuaternion(quaternion);
        }

        const theta = Math.acos(Math.max(-1, Math.min(1, physDir.z)));
        const phi = Math.atan2(physDir.y, physDir.x);

        let rMin = 0.1;
        let rMax = baseRadius * 1.5;
        let radius95 = rMax * 0.5;

        const densityAtMax = calcDensity(rMax, theta, phi);
        if (densityAtMax > densityThreshold) {
            rMax = baseRadius * 2.5;
        }

        for (let iter = 0; iter < 25; iter++) {
            const rMid = (rMin + rMax) / 2;
            const density = calcDensity(rMid, theta, phi);

            if (density > densityThreshold) {
                rMin = rMid;
            } else {
                rMax = rMid;
            }
            radius95 = rMid;
        }

        vertices[i] = physDir.x * radius95;
        vertices[i + 1] = physDir.y * radius95;
        vertices[i + 2] = physDir.z * radius95;

        colors[i] = 0.2;
        colors[i + 1] = 1.0;
        colors[i + 2] = 0.5;
    }

    geometry.setAttribute('color', new THREE.BufferAttribute(colors, 3));
    geometry.attributes.position.needsUpdate = true;

    const material = new THREE.MeshBasicMaterial({
        vertexColors: true,
        transparent: true,
        opacity: 0.6,
        side: THREE.DoubleSide,
        wireframe: true,
        wireframeLinewidth: 1.5
    });

    const mesh = new THREE.Mesh(geometry, material);
    mesh.layers.set(1);
    group.add(mesh);
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
 * 只为可见的轨道创建等值面
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

    // 对轨道进行排序
    if (Hydrogen.sortOrbitalsForHybridization) {
        orbitalParams = Hydrogen.sortOrbitalsForHybridization(orbitalParams);
    }

    const numOrbitals = orbitalParams.length;
    const coeffMatrix = Hydrogen.getHybridCoefficients ? Hydrogen.getHybridCoefficients(numOrbitals) : null;

    if (!coeffMatrix) return overlays;

    // 不同颜色用于区分各个杂化轨道
    const colors = [
        { r: 0.2, g: 1.0, b: 0.5 },   // 绿色
        { r: 1.0, g: 0.5, b: 0.2 },   // 橙色
        { r: 0.2, g: 0.5, b: 1.0 },   // 蓝色
        { r: 1.0, g: 0.2, b: 0.8 },   // 粉色
        { r: 0.8, g: 1.0, b: 0.2 },   // 黄绿色
        { r: 0.2, g: 1.0, b: 1.0 },   // 青色
    ];

    // 为每个可见的杂化轨道创建等值面
    for (let hybridIndex = 0; hybridIndex < numOrbitals; hybridIndex++) {
        // 检查该轨道是否可见
        if (state.hybridOrbitalVisibility && state.hybridOrbitalVisibility[hybridIndex] === false) {
            overlays.push(null); // 占位，保持索引一致
            continue;
        }

        const pointIndices = state.hybridOrbitalPointsMap[hybridIndex] || [];
        if (pointIndices.length < 50) {
            overlays.push(null);
            continue;
        }

        const overlay = createSingleHybridContour(
            orbitalParams,
            coeffMatrix[hybridIndex],
            hybridIndex,
            pointIndices,
            colors[hybridIndex % colors.length]
        );

        if (overlay) {
            overlays.push(overlay);
        } else {
            overlays.push(null);
        }
    }

    return overlays;

    // 内部函数：创建单个杂化轨道的等值面
    function createSingleHybridContour(orbitalParams, coeffs, hybridIndex, pointIndices, color) {
        const overlayGroup = new THREE.Group();

        // 获取该杂化轨道的主轴方向
        let symmetryAxis = { x: 0, y: 0, z: 1 };
        if (Hydrogen.getHybridLobeAxis) {
            symmetryAxis = Hydrogen.getHybridLobeAxis(orbitalParams, hybridIndex);
        } else if (Hydrogen.findHybridPrincipalAxis) {
            symmetryAxis = Hydrogen.findHybridPrincipalAxis(orbitalParams, hybridIndex);
        }

        // 计算旋转四元数
        let quaternion = null;
        const threeJsPoleAxis = new THREE.Vector3(0, 1, 0);
        const targetAxis = new THREE.Vector3(symmetryAxis.x, symmetryAxis.y, symmetryAxis.z).normalize();

        const dotProduct = threeJsPoleAxis.dot(targetAxis);
        if (Math.abs(dotProduct - 1) > 1e-6) {
            quaternion = new THREE.Quaternion();
            quaternion.setFromUnitVectors(threeJsPoleAxis, targetAxis);
        }

        // 密度计算函数
        function calcDensity(r, theta, phi) {
            let psi = 0;
            for (let j = 0; j < orbitalParams.length; j++) {
                const op = orbitalParams[j];
                const R = Hydrogen.radialWavefunction(op.n, op.l, r);
                const Y = Hydrogen.realYlm_value(op.angKey.l, op.angKey.m, op.angKey.t, theta, phi);
                psi += coeffs[j] * R * Y;
            }
            return psi * psi;
        }

        // 计算该轨道点的95%分位密度阈值
        const positions = state.points.geometry.attributes.position.array;
        const densities = [];

        for (const pointIdx of pointIndices) {
            const i3 = pointIdx * 3;
            const x = positions[i3];
            const y = positions[i3 + 1];
            const z = positions[i3 + 2];
            const r = Math.sqrt(x * x + y * y + z * z);
            if (r < 1e-10) continue;

            const theta = Math.acos(Math.max(-1, Math.min(1, z / r)));
            const phi = Math.atan2(y, x);
            const density = calcDensity(r, theta, phi);
            densities.push(density);
        }

        if (densities.length < 50) return null;

        densities.sort((a, b) => b - a);
        const index95 = Math.floor(densities.length * 0.95);
        const densityThreshold = densities[index95] || 0;

        // 估算平均半径
        let sumR = 0;
        for (const pointIdx of pointIndices) {
            const i3 = pointIdx * 3;
            const x = positions[i3];
            const y = positions[i3 + 1];
            const z = positions[i3 + 2];
            sumR += Math.sqrt(x * x + y * y + z * z);
        }
        const avgRadius = sumR / pointIndices.length;

        // 创建网格
        const thetaBins = 32;
        const phiBins = 64;
        const geometry = new THREE.SphereGeometry(1, phiBins, thetaBins);
        const vertices = geometry.attributes.position.array;
        const colorsAttr = new Float32Array(vertices.length);

        for (let i = 0; i < vertices.length; i += 3) {
            const x0 = vertices[i];
            const y0 = vertices[i + 1];
            const z0 = vertices[i + 2];

            const unitR = Math.sqrt(x0 * x0 + y0 * y0 + z0 * z0);
            if (unitR < 1e-10) continue;

            let physDir = new THREE.Vector3(x0 / unitR, y0 / unitR, z0 / unitR);
            if (quaternion) {
                physDir.applyQuaternion(quaternion);
            }

            const theta = Math.acos(Math.max(-1, Math.min(1, physDir.z)));
            const phi = Math.atan2(physDir.y, physDir.x);

            // 二分搜索找到95%等值面半径
            let rMin = 0.1;
            let rMax = avgRadius * 2.0;
            let radius95 = avgRadius;

            for (let iter = 0; iter < 20; iter++) {
                const rMid = (rMin + rMax) / 2;
                const density = calcDensity(rMid, theta, phi);

                if (density > densityThreshold) {
                    rMin = rMid;
                } else {
                    rMax = rMid;
                }
                radius95 = rMid;
            }

            vertices[i] = physDir.x * radius95;
            vertices[i + 1] = physDir.y * radius95;
            vertices[i + 2] = physDir.z * radius95;

            colorsAttr[i] = color.r;
            colorsAttr[i + 1] = color.g;
            colorsAttr[i + 2] = color.b;
        }

        geometry.setAttribute('color', new THREE.BufferAttribute(colorsAttr, 3));
        geometry.attributes.position.needsUpdate = true;

        const material = new THREE.MeshBasicMaterial({
            vertexColors: true,
            transparent: true,
            opacity: 0.5,
            side: THREE.DoubleSide,
            wireframe: true,
            wireframeLinewidth: 1.5
        });

        const mesh = new THREE.Mesh(geometry, material);
        mesh.layers.set(1);
        overlayGroup.add(mesh);

        return overlayGroup;
    }
};
