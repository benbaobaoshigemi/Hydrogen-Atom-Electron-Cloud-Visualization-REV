// UI 控件交互模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.UI = {};

// 防抖函数，避免频繁触发
window.ElectronCloud.UI.debounce = function (func, wait) {
    let timeout;
    return function executedFunction(...args) {
        const later = () => {
            clearTimeout(timeout);
            func(...args);
        };
        clearTimeout(timeout);
        timeout = setTimeout(later, wait);
    };
};

// 初始化UI事件监听
window.ElectronCloud.UI.init = function () {
    const ui = window.ElectronCloud.ui;

    // 初始化顶部模式切换栏（必须在最开始调用，绑定点击事件）
    if (window.ElectronCloud.UI.initModeSwitcher) {
        window.ElectronCloud.UI.initModeSwitcher();
    }

    // 初始化杂化轨道控制面板
    if (window.ElectronCloud.UI.initHybridPanel) {
        window.ElectronCloud.UI.initHybridPanel();
    }

    // 初始化自定义轨道列表
    if (window.ElectronCloud.UI.initCustomOrbitalList) {
        window.ElectronCloud.UI.initCustomOrbitalList();
    }

    // 动态填充原子列表（必须在 initCustomSelects 之前调用，否则自定义UI会包含旧选项）
    // populateAtomList(); // DISABLED: Using hardcoded list in HTML now

    // 初始化通用自定义下拉框
    if (window.ElectronCloud.UI.initCustomSelects) {
        window.ElectronCloud.UI.initCustomSelects();
    }

    // 监听UI事件
    if (ui.axesSizeRange) {
        ui.axesSizeRange.addEventListener('input', window.ElectronCloud.UI.onAxesSizeChange);
    }

    if (ui.angular3dToggle) {
        ui.angular3dToggle.addEventListener('change', window.ElectronCloud.UI.onAngularOverlayToggle);
    }

    // 95%等值面轮廓开关
    const contour3dToggle = document.getElementById('contour-3d-toggle');
    if (contour3dToggle) {
        contour3dToggle.addEventListener('change', window.ElectronCloud.UI.onContourOverlayToggle);
    }

    // 杂化模式轨道轮廓开关
    const hybridContourToggle = document.getElementById('hybrid-contour-toggle');
    if (hybridContourToggle) {
        hybridContourToggle.addEventListener('change', window.ElectronCloud.UI.onHybridContourToggle);
    }

    if (ui.orbitalSelect) {
        ui.orbitalSelect.addEventListener('change', window.ElectronCloud.UI.onOrbitalSelectChange);
    }

    if (ui.startButton) {
        ui.startButton.addEventListener('click', window.ElectronCloud.Orbital.startDrawing);
    }

    if (ui.speedRange) {
        ui.speedRange.addEventListener('input', window.ElectronCloud.UI.onSpeedChange);
    }

    if (ui.sizeRange) {
        ui.sizeRange.addEventListener('input', window.ElectronCloud.UI.onSizeChange);
    }

    if (ui.opacityRange) {
        ui.opacityRange.addEventListener('input', window.ElectronCloud.UI.onOpacityChange);
    }

    if (ui.centerLock) {
        ui.centerLock.addEventListener('change', window.ElectronCloud.UI.onCenterLockChange);
    }

    if (ui.pauseButton) {
        ui.pauseButton.addEventListener('click', window.ElectronCloud.UI.onPauseToggle);
    }

    if (ui.clearButton) {
        ui.clearButton.addEventListener('click', () => {
            const state = window.ElectronCloud.state;

            // 【用户需求】在多选模式和杂化模式下，重置不应清除轨道选择
            // 保存当前模式和选择
            const savedOrbitals = [...state.currentOrbitals];
            const savedIsHybridMode = state.isHybridMode;
            const isMultiOrHybrid = (savedOrbitals.length > 1 || savedIsHybridMode);

            // 重置状态（清除采样数据等）
            window.ElectronCloud.resetState();
            window.ElectronCloud.Orbital.clearDrawing();

            // 如果是多选或杂化模式，恢复轨道选择
            if (isMultiOrHybrid && savedOrbitals.length > 0) {
                state.currentOrbitals = savedOrbitals;
                state.isHybridMode = savedIsHybridMode;
                // 恢复UI状态
                if (ui.orbitalSelect) {
                    Array.from(ui.orbitalSelect.options).forEach(option => {
                        if (savedOrbitals.includes(option.value)) {
                            option.classList.add('force-selected');
                            // 【关键修复】确保原生 select 的选中状态也同步，
                            // 这样 updateHybridOrbitalButtons 才能正确读取 selectedOptions
                            option.selected = true;
                        }
                    });
                }
                window.ElectronCloud.UI.updateCustomListVisuals();
                window.ElectronCloud.UI.updateSelectionCount();

                // 【关键修复】如果恢复了杂化模式，需要重置并刷新杂化面板（生成按钮等）
                if (savedIsHybridMode && window.ElectronCloud.UI.resetHybridPanel) {
                    window.ElectronCloud.UI.resetHybridPanel();
                }
            } else {
                // 【新增】单选模式下，重置时清除所有选择
                // 确保用户需要重新选择轨道才能启动渲染
                if (ui.orbitalSelect) {
                    Array.from(ui.orbitalSelect.options).forEach(option => {
                        option.selected = false;
                    });
                }
                state.currentOrbitals = [];
                window.ElectronCloud.UI.updateCustomListVisuals();
            }

            if (window.DataPanel && window.DataPanel.reset) {
                window.DataPanel.reset();
            }
            // 更新模式切换栏状态
            window.ElectronCloud.UI.updateModeSwitcherState();
            // 更新3D角向形状和轮廓开关状态（因pointCount归零，应自动禁用）
            if (window.ElectronCloud.UI.updateAngular3DToggleState) {
                window.ElectronCloud.UI.updateAngular3DToggleState();
            }
        });
    }

    // 实验功能面板折叠逻辑
    const experimentalPanel = document.getElementById('experimental-panel');
    const experimentalCollapseBtn = document.getElementById('experimental-collapse-btn');
    const experimentalPanelTab = document.getElementById('experimental-panel-tab');

    if (experimentalPanel && experimentalCollapseBtn && experimentalPanelTab) {
        experimentalCollapseBtn.addEventListener('click', () => {
            experimentalPanel.classList.add('collapsed');
        });

        experimentalPanelTab.addEventListener('click', () => {
            experimentalPanel.classList.remove('collapsed');
        });
    }

    // 手势状态图标点击逻辑 (停止手势)
    const gestureStatusIcon = document.getElementById('gesture-status-icon');
    if (gestureStatusIcon) {
        gestureStatusIcon.addEventListener('click', () => {
            if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.stop) {
                window.ElectronCloud.Gesture.stop();
                // 退出全屏
                if (document.fullscreenElement) {
                    document.exitFullscreen().catch(err => console.error(err));
                }
            }
        });
    }

    // 全屏按钮逻辑
    const fullscreenBtn = document.getElementById('app-fullscreen-btn');
    if (fullscreenBtn) {
        fullscreenBtn.addEventListener('click', () => {
            if (!document.fullscreenElement) {
                document.documentElement.requestFullscreen().catch(err => {
                    console.error(`Error attempting to enable full-screen mode: ${err.message} (${err.name})`);
                });
            } else {
                document.exitFullscreen();
            }
        });

        // 监听全屏状态变化更新图标
        document.addEventListener('fullscreenchange', () => {
            if (document.fullscreenElement) {
                fullscreenBtn.textContent = '⛶'; // 或者换成退出全屏的图标
                fullscreenBtn.title = '退出全屏';
                fullscreenBtn.classList.add('active');
            } else {
                fullscreenBtn.textContent = '⛶';
                fullscreenBtn.title = '全屏模式';
                fullscreenBtn.classList.remove('active');
            }
        });
    }

    // 手势控制按钮逻辑
    const gestureBtn = document.getElementById('gesture-control-btn');
    if (gestureBtn) {
        gestureBtn.addEventListener('click', window.ElectronCloud.UI.onGestureButtonClick);
    }

    if (ui.plotTypeSelect) {
        ui.plotTypeSelect.addEventListener('change', () => {
            window.ElectronCloud.Orbital.drawProbabilityChart();
        });
    }

    // 多选模式功能
    if (ui.multiselectToggle && ui.orbitalSelect) {
        ui.multiselectToggle.addEventListener('change', window.ElectronCloud.UI.onMultiselectToggle);
    }

    // 比照模式功能
    if (ui.compareToggle && ui.orbitalSelect) {
        ui.compareToggle.addEventListener('change', window.ElectronCloud.UI.onCompareToggle);
    }

    // 设置多选模式的交互
    if (ui.orbitalSelect) {
        window.ElectronCloud.UI.setupMultiselectMode();
    }

    // 清除所有选择按钮
    if (ui.clearAllSelectionsBtn) {
        ui.clearAllSelectionsBtn.addEventListener('click', window.ElectronCloud.UI.clearAllSelections);
    }

    // 为所有 mode-toggle-box 添加点击切换逻辑
    window.ElectronCloud.UI.setupModeToggleBoxes();

    // 强制同步一次多选状态
    window.ElectronCloud.UI.onMultiselectToggle();

    // 初始化单选模式的选中样式
    window.ElectronCloud.UI.updateSingleSelectStyle();

    // 确保初始状态下显示开关的状态正确
    window.ElectronCloud.UI.updateAngular3DToggleState();

    // 确保初始状态下清空按钮状态正确
    window.ElectronCloud.UI.updateClearAllSelectionsState();

    // 确保初始状态下坐标轴滑动条状态正确
    window.ElectronCloud.UI.updateAxesSizeRangeState();

    // 确保初始状态下模式切换栏状态正确
    window.ElectronCloud.UI.updateModeSwitcherState();

    // 确保初始状态下自动旋转按钮可用（允许预设）
    window.ElectronCloud.UI.updateAutoRotateButtonState();

    // 确保初始状态下坐标系正确隐藏（farthestDistance=0时）
    if (window.ElectronCloud.ui.axesSizeRange) {
        window.ElectronCloud.UI.onAxesSizeChange({ target: window.ElectronCloud.ui.axesSizeRange });
    }

    // 【性能优化】比照模式选择器采用懒加载，首次切换到比照模式时才初始化
    // 避免启动时同步创建312个自定义DOM节点阻塞主线程
    // window.ElectronCloud.UI.initCompareOrbitalSelectors(); // LAZY LOAD
    window.ElectronCloud.UI._compareSelectorsInitialized = false;

    // 【性能优化】使用事件委托统一处理所有自定义下拉框的关闭
    document.addEventListener('click', window.ElectronCloud.UI._handleGlobalClickForCustomSelects);

    // 视图面板折叠逻辑
    const viewPanel = document.getElementById('view-panel');
    const viewCollapseBtn = document.getElementById('view-collapse-btn');
    const viewPanelTab = document.getElementById('view-panel-tab');

    if (viewPanel && viewCollapseBtn && viewPanelTab) {
        viewCollapseBtn.addEventListener('click', () => {
            viewPanel.classList.add('collapsed');
        });

        viewPanelTab.addEventListener('click', () => {
            viewPanel.classList.remove('collapsed');
        });
    }

    // 旋转轴输入
    const rotationAxisX = document.getElementById('rotation-axis-x');
    const rotationAxisY = document.getElementById('rotation-axis-y');
    const rotationAxisZ = document.getElementById('rotation-axis-z');
    const rotationSpeedRange = document.getElementById('rotation-speed-range');
    const rotationSpeedValue = document.getElementById('rotation-speed-value');

    const updateRotationAxis = () => {
        const state = window.ElectronCloud.state;
        const x = parseFloat(rotationAxisX?.value) || 0;
        const y = parseFloat(rotationAxisY?.value) || 0;
        const z = parseFloat(rotationAxisZ?.value) || 0;

        state.autoRotate.axis.set(x, y, z);

        // 如果有有效轴，显示辅助线
        if (x !== 0 || y !== 0 || z !== 0) {
            if (window.ElectronCloud.Scene.showRotationAxisHelper) {
                window.ElectronCloud.Scene.showRotationAxisHelper(state.autoRotate.axis);
            }
        }
    };

    if (rotationAxisX) rotationAxisX.addEventListener('change', updateRotationAxis);
    if (rotationAxisY) rotationAxisY.addEventListener('change', updateRotationAxis);
    if (rotationAxisZ) rotationAxisZ.addEventListener('change', updateRotationAxis);

    // 旋转速度滑动条（只设置速度，不自动启动）
    if (rotationSpeedRange) {
        rotationSpeedRange.addEventListener('input', (e) => {
            const state = window.ElectronCloud.state;
            const value = parseFloat(e.target.value);
            state.autoRotate.speed = value;
            // 不再自动启动，由开关控制

            if (rotationSpeedValue) {
                rotationSpeedValue.textContent = value.toFixed(2);
            }
        });
    }

    // 自动旋转开关
    const autoRotateToggle = document.getElementById('auto-rotate-toggle');
    const rotationFeatureBox = document.getElementById('rotation-feature-box');
    const recordRotationBtn = document.getElementById('record-rotation-btn');

    if (autoRotateToggle) {
        autoRotateToggle.addEventListener('change', (e) => {
            const state = window.ElectronCloud.state;
            state.autoRotate.enabled = e.target.checked;

            // 功能框激活状态变绿
            if (rotationFeatureBox) {
                rotationFeatureBox.classList.toggle('active', e.target.checked);
            }

            // 更新录制按钮状态
            window.ElectronCloud.UI.updateRecordButtonState();

            // 互斥逻辑：自动旋转启用时，禁用锁定视角
            window.ElectronCloud.UI.updateRotationLockMutualExclusion();
        });
    }

    // 点击 rotation-feature-box 切换开关
    if (rotationFeatureBox) {
        rotationFeatureBox.addEventListener('click', (e) => {
            if (autoRotateToggle && !autoRotateToggle.disabled) {
                autoRotateToggle.checked = !autoRotateToggle.checked;
                autoRotateToggle.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    }

    // 录制旋转视频按钮
    if (recordRotationBtn) {
        recordRotationBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            window.ElectronCloud.UI.toggleRotationRecording();
        });
    }

    // 闪烁模式开关
    const heartbeatToggle = document.getElementById('heartbeat-toggle');
    const flickerSettingsBtn = document.getElementById('flicker-settings-btn');
    const flickerSubmenu = document.getElementById('flicker-submenu');
    const flickerFeatureBox = document.getElementById('flicker-feature-box');
    const flickerModeSelect = document.getElementById('flicker-mode-select');
    const flickerFrequencyRange = document.getElementById('flicker-frequency-range');
    const heartbeatMaxBrightness = document.getElementById('heartbeat-max-brightness');

    if (heartbeatToggle) {
        heartbeatToggle.addEventListener('change', (e) => {
            const state = window.ElectronCloud.state;
            state.heartbeat.enabled = e.target.checked;

            // 功能框激活状态变绿
            if (flickerFeatureBox) {
                flickerFeatureBox.classList.toggle('active', e.target.checked);
            }

            // 闪烁模式下仍然允许调节亮度（不再禁用滑动条）
            // 亮度滑动条控制点云基础透明度，闪烁效果在此基础上叠加

            // 如果关闭闪烁模式，恢复bloom状态、原始颜色和透明度
            if (!e.target.checked && state.bloomPass) {
                // 【修复】根据亮度滑块当前值恢复辉光状态
                const ui = window.ElectronCloud.ui;
                const currentValue = ui.opacityRange ? parseInt(ui.opacityRange.value, 10) : 67;
                if (currentValue > 50) {
                    // 辉光区间：50-100 映射到辉光强度 0-4.5
                    const glowProgress = (currentValue - 50) / 50;
                    const curvedProgress = glowProgress * glowProgress;
                    state.bloomPass.strength = curvedProgress * 4.5;
                    state.bloomEnabled = true;
                } else {
                    // 透明度区间，关闭辉光
                    state.bloomPass.strength = 0;
                    state.bloomEnabled = false;
                }

                // 【修复】根据亮度滑块当前值恢复透明度
                if (state.points && state.points.material) {
                    if (currentValue <= 50) {
                        // 透明度区间：0-50 映射到 0.05-1.0
                        state.points.material.opacity = 0.05 + (currentValue / 50) * 0.95;
                    } else {
                        // 辉光区间：透明度保持最大
                        state.points.material.opacity = 1.0;
                    }
                    state.points.material.needsUpdate = true;
                }

                // 恢复星空模式下修改的原始颜色
                if (state.baseColors && state.points && state.points.geometry) {
                    const colors = state.points.geometry.attributes.color;
                    if (colors) {
                        const baseColors = state.baseColors;
                        const pointCount = state.pointCount || 0;
                        for (let i = 0; i < pointCount * 3 && i < baseColors.length; i++) {
                            colors.array[i] = baseColors[i];
                        }
                        colors.needsUpdate = true;
                    }
                }
                // 清除星空模式的缓存数据（保留baseColors以便重新开启时使用）
                state.heartbeat.starryPhases = null;
                state.heartbeat.starryFrequencies = null;
            }
        });
    } else {
        console.error('heartbeat-toggle 元素未找到');
    }

    // 点击 flicker-feature-box 切换闪烁模式
    if (flickerFeatureBox) {
        flickerFeatureBox.addEventListener('click', (e) => {
            if (heartbeatToggle) {
                heartbeatToggle.checked = !heartbeatToggle.checked;
                heartbeatToggle.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    }

    // 闪烁设置二级面板
    if (flickerSettingsBtn && flickerSubmenu) {
        flickerSettingsBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            const isVisible = flickerSubmenu.classList.contains('visible');
            closeAllSubmenus();
            if (!isVisible) {
                // 获取实验面板的位置
                const experimentalPanel = document.getElementById('experimental-panel');
                const panelRect = experimentalPanel ? experimentalPanel.getBoundingClientRect() : null;
                const btnRect = flickerSettingsBtn.getBoundingClientRect();

                // 二级菜单定位在面板左侧，完全不重叠
                flickerSubmenu.style.top = btnRect.top + 'px';
                flickerSubmenu.style.right = 'auto';
                flickerSubmenu.style.left = 'auto';
                if (panelRect) {
                    // 定位在面板左边缘的左侧，留出10px间距
                    flickerSubmenu.style.right = (window.innerWidth - panelRect.left + 10) + 'px';
                } else {
                    flickerSubmenu.style.right = (window.innerWidth - btnRect.left + 10) + 'px';
                }
                flickerSubmenu.classList.add('visible');
                flickerSettingsBtn.classList.add('active');
            }
        });
    }

    // 闪烁模式选择
    if (flickerModeSelect) {
        flickerModeSelect.addEventListener('change', (e) => {
            const state = window.ElectronCloud.state;
            const previousMode = state.heartbeat.mode;
            state.heartbeat.mode = e.target.value;
            state.heartbeat.phase = 0;

            // 【重要】切换模式时，先恢复原始颜色，避免颜色污染
            if (state.baseColors && state.points && state.points.geometry) {
                const colors = state.points.geometry.attributes.color;
                if (colors) {
                    const pointCount = state.pointCount || 0;
                    for (let i = 0; i < pointCount * 3 && i < state.baseColors.length; i++) {
                        colors.array[i] = state.baseColors[i];
                    }
                    colors.needsUpdate = true;
                }
            }

            // 如果从星空模式切换出去，清除星空缓存
            if (previousMode === 'starry') {
                state.heartbeat.starryPhases = null;
                state.heartbeat.starryFrequencies = null;
                state.heartbeat.originalColors = null;
            }

            // 切换到新模式时，重置该模式的缓存（强制重新计算）
            // 这解决了"第一次启用是旧样式"的问题
            state.diffuseDensitiesComputed = false;
            state.waveRanksComputed = false;
            state.waveRanksPointCount = 0;
        });
    }

    // 闪烁频率
    if (flickerFrequencyRange) {
        flickerFrequencyRange.addEventListener('input', (e) => {
            const state = window.ElectronCloud.state;
            state.heartbeat.frequency = parseInt(e.target.value, 10);
        });
    }

    // 最大亮度滑动条
    if (heartbeatMaxBrightness) {
        heartbeatMaxBrightness.addEventListener('input', (e) => {
            const state = window.ElectronCloud.state;
            state.heartbeat.maxBrightness = parseInt(e.target.value, 10);
        });
    }

    // 权重模式开关 - 独立实验功能
    const weightModeToggle = document.getElementById('weight-mode-toggle');
    const weightFeatureBox = document.getElementById('weight-feature-box');

    if (weightModeToggle) {
        weightModeToggle.addEventListener('change', (e) => {
            const state = window.ElectronCloud.state;
            state.weightMode = e.target.checked;
            // 兼容旧的heartbeat.weightMode
            if (state.heartbeat) {
                state.heartbeat.weightMode = e.target.checked;
            }

            // 更新框的激活状态
            const box = document.getElementById('weight-feature-box');
            if (box) {
                box.classList.toggle('active', e.target.checked);
            }

            // 切换权重模式时，清除缓存的密度权重，以便重新计算
            if (!e.target.checked) {
                state.densityWeights = null;
            }
        });
    } else {
        console.error('weight-mode-toggle 元素未找到');
    }

    // 权重模式功能框点击切换
    if (weightFeatureBox) {
        weightFeatureBox.addEventListener('click', (e) => {
            const toggle = document.getElementById('weight-mode-toggle');
            if (toggle) {
                toggle.checked = !toggle.checked;
                toggle.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    } else {
        console.error('weight-feature-box 元素未找到');
    }

    // 滚动生成模式开关
    const rollingModeToggle = document.getElementById('rolling-mode-toggle');
    const rollingFeatureBox = document.getElementById('rolling-feature-box');

    if (rollingModeToggle) {
        rollingModeToggle.addEventListener('change', (e) => {
            const state = window.ElectronCloud.state;

            // 渲染完成后直接生效；渲染前只设置pending
            if (state.samplingCompleted) {
                state.rollingMode.enabled = e.target.checked;
            } else {
                state.rollingMode.pendingEnabled = e.target.checked;
            }

            // 更新按钮激活样式
            if (rollingFeatureBox) {
                rollingFeatureBox.classList.toggle('active', e.target.checked);
            }
        });
    }

    if (rollingFeatureBox) {
        rollingFeatureBox.addEventListener('click', (e) => {
            if (rollingModeToggle) {
                rollingModeToggle.checked = !rollingModeToggle.checked;
                rollingModeToggle.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    }

    // 导出设置二级面板
    const exportSettingsBtn = document.getElementById('export-settings-btn');
    const exportSubmenu = document.getElementById('export-submenu');
    const exportFeatureBox = document.getElementById('export-feature-box');

    if (exportSettingsBtn && exportSubmenu) {
        exportSettingsBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            const isVisible = exportSubmenu.classList.contains('visible');
            closeAllSubmenus();
            if (!isVisible) {
                // 获取实验面板的位置
                const experimentalPanel = document.getElementById('experimental-panel');
                const panelRect = experimentalPanel ? experimentalPanel.getBoundingClientRect() : null;
                const btnRect = exportSettingsBtn.getBoundingClientRect();

                // 二级菜单定位在面板左侧，完全不重叠
                exportSubmenu.style.top = btnRect.top + 'px';
                exportSubmenu.style.right = 'auto';
                exportSubmenu.style.left = 'auto';
                if (panelRect) {
                    // 定位在面板左边缘的左侧，留出10px间距
                    exportSubmenu.style.right = (window.innerWidth - panelRect.left + 10) + 'px';
                } else {
                    exportSubmenu.style.right = (window.innerWidth - btnRect.left + 10) + 'px';
                }
                exportSubmenu.classList.add('visible');
                exportSettingsBtn.classList.add('active');
            }
        });
    }

    // 点击导出功能框直接导出图片
    if (exportFeatureBox) {
        exportFeatureBox.addEventListener('click', (e) => {
            if (window.ElectronCloud.Visualization && typeof window.ElectronCloud.Visualization.exportImage === 'function') {
                window.ElectronCloud.Visualization.exportImage();
            } else {
                console.error('导出函数不可用');
                alert('导出功能暂时不可用，请刷新页面后重试');
            }
        });
    }

    // 旋转设置二级面板
    const rotationSettingsBtn = document.getElementById('rotation-settings-btn');
    const rotationSubmenu = document.getElementById('rotation-submenu');

    if (rotationSettingsBtn && rotationSubmenu) {
        rotationSettingsBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            const isVisible = rotationSubmenu.classList.contains('visible');
            closeAllSubmenus();
            if (!isVisible) {
                // 获取视图面板的位置
                const viewPanel = document.getElementById('view-panel');
                const panelRect = viewPanel ? viewPanel.getBoundingClientRect() : null;
                const btnRect = rotationSettingsBtn.getBoundingClientRect();

                // 二级菜单定位在面板右侧，完全不重叠
                rotationSubmenu.style.top = btnRect.top + 'px';
                rotationSubmenu.style.right = 'auto';
                rotationSubmenu.style.left = 'auto';
                if (panelRect) {
                    // 定位在面板右边缘的右侧，留出10px间距
                    rotationSubmenu.style.left = (panelRect.right + 10) + 'px';
                } else {
                    rotationSubmenu.style.left = (btnRect.right + 10) + 'px';
                }
                rotationSubmenu.classList.add('visible');
                rotationSettingsBtn.classList.add('active');
            }
        });
    }

    // 视角锁定按钮
    const viewLockBtns = document.querySelectorAll('.view-lock-btn');
    viewLockBtns.forEach(btn => {
        btn.addEventListener('click', () => {
            const axis = btn.dataset.axis;
            const state = window.ElectronCloud.state;

            // 切换激活状态
            const wasActive = btn.classList.contains('active');
            // 清除 active 和 green 状态
            viewLockBtns.forEach(b => { b.classList.remove('active'); b.classList.remove('green'); });

            if (!wasActive) {
                btn.classList.add('active');
                btn.classList.add('green');
                // 锁定到对应视角
                window.ElectronCloud.UI.lockCameraToAxis(axis);
                // 互斥逻辑：锁定视角时禁用自动旋转
                window.ElectronCloud.UI.updateRotationLockMutualExclusion();
            } else {
                // 解锁 - 使用统一的旋转重置函数
                // 这是解决"坐标歪斜"问题的关键：确保解锁时所有对象旋转同步归零
                if (window.ElectronCloud.Scene.resetAllSceneObjectsRotation) {
                    window.ElectronCloud.Scene.resetAllSceneObjectsRotation();
                }

                // 解除锁定时更新互斥状态
                window.ElectronCloud.UI.updateRotationLockMutualExclusion();
            }
        });
    });

    // 关闭所有子菜单的函数
    function closeAllSubmenus() {
        document.querySelectorAll('.submenu-panel').forEach(panel => {
            panel.classList.remove('visible');
        });
        document.querySelectorAll('.feature-settings-btn, .submenu-trigger-right').forEach(btn => {
            btn.classList.remove('active');
        });
    }

    // 点击其他地方关闭子菜单
    document.addEventListener('click', (e) => {
        if (!e.target.closest('.submenu-panel') &&
            !e.target.closest('.feature-settings-btn') &&
            !e.target.closest('.submenu-trigger-right')) {
            closeAllSubmenus();
        }
    });

    // 初始化亮度/辉光状态（基于默认滑块值）
    // 延迟执行，确保scene.js已初始化bloomPass
    setTimeout(() => {
        if (window.ElectronCloud.UI.onOpacityChange) {
            window.ElectronCloud.UI.onOpacityChange();
        }
    }, 100);

    // 相位显示开关
    const phaseToggle = document.getElementById('phase-toggle');
    if (phaseToggle) {
        phaseToggle.addEventListener('change', window.ElectronCloud.UI.onPhaseToggle);
    }

    // ==================== 原子选择器 ====================
    const atomSettingsBtn = document.getElementById('atom-settings-btn');
    const atomSubmenu = document.getElementById('atom-submenu');
    const atomFeatureBox = document.getElementById('atom-feature-box');

    if (atomSettingsBtn && atomSubmenu) {
        // 设置按钮点击：展开/收起子面板
        atomSettingsBtn.addEventListener('click', (e) => {
            e.stopPropagation();
            const state = window.ElectronCloud.state;

            // 渲染中或已完成时禁用
            if (state.isDrawing || state.pointCount > 0) {
                return;
            }

            // 切换子面板显示
            const isVisible = atomSubmenu.classList.contains('visible');
            closeAllSubmenus();
            if (!isVisible) {
                // 获取实验面板的位置（与闪烁模式相同的定位逻辑）
                const experimentalPanel = document.getElementById('experimental-panel');
                const panelRect = experimentalPanel ? experimentalPanel.getBoundingClientRect() : null;
                const btnRect = atomSettingsBtn.getBoundingClientRect();

                // 二级菜单定位在面板左侧，完全不重叠
                atomSubmenu.style.top = btnRect.top + 'px';
                atomSubmenu.style.right = 'auto';
                atomSubmenu.style.left = 'auto';
                if (panelRect) {
                    // 定位在面板左边缘的左侧，留出10px间距
                    atomSubmenu.style.right = (window.innerWidth - panelRect.left + 10) + 'px';
                } else {
                    atomSubmenu.style.right = (window.innerWidth - btnRect.left + 10) + 'px';
                }
                atomSubmenu.classList.add('visible');
                atomSettingsBtn.classList.add('active');
            }
        });
    }

    // 原子下拉选择框事件
    const atomSelect = document.getElementById('atom-select');
    if (atomSelect) {
        atomSelect.addEventListener('change', function () {
            const atomSymbol = this.value;
            const state = window.ElectronCloud.state;

            // 渲染中或已完成时禁用
            if (state.isDrawing || state.pointCount > 0) {
                this.value = state.currentAtom; // 恢复原值
                return;
            }

            // 更新显示
            const currentLabel = document.getElementById('current-atom-label');
            if (currentLabel) currentLabel.textContent = atomSymbol;

            // 更新状态
            state.currentAtom = atomSymbol;

            // 更新原子信息显示
            const atomInfo = document.getElementById('atom-info');
            let groundState = '1s¹'; // 默认H

            if (atomSymbol === 'H') {
                groundState = '1s¹';
            } else if (window.SlaterBasis && window.SlaterBasis[atomSymbol]) {
                const atom = window.SlaterBasis[atomSymbol];
                groundState = atom.groundState || '';
            }
            if (atomInfo) atomInfo.textContent = `基态: ${groundState}`;

            // 更新轨道选择器为该原子可用的轨道
            if (window.ElectronCloud.UI.updateOrbitalSelector) {
                window.ElectronCloud.UI.updateOrbitalSelector(atomSymbol);
            }
            console.log('选择原子:', atomSymbol);
        });
    }

    // 动态填充原子选择列表
    function populateAtomList() {
        const select = document.getElementById('atom-select');
        if (!select || !window.SlaterBasis) return;

        // 元素中文名称映射 (1-54)
        const cnNames = {
            'H': '氢', 'He': '氦', 'Li': '锂', 'Be': '铍', 'B': '硼', 'C': '碳', 'N': '氮', 'O': '氧', 'F': '氟', 'Ne': '氖',
            'Na': '钠', 'Mg': '镁', 'Al': '铝', 'Si': '硅', 'P': '磷', 'S': '硫', 'Cl': '氯', 'Ar': '氩',
            'K': '钾', 'Ca': '钙', 'Sc': '钪', 'Ti': '钛', 'V': '钒', 'Cr': '铬', 'Mn': '锰', 'Fe': '铁', 'Co': '钴', 'Ni': '镍', 'Cu': '铜', 'Zn': '锌',
            'Ga': '镓', 'Ge': '锗', 'As': '砷', 'Se': '硒', 'Br': '溴', 'Kr': '氪',
            'Rb': '铷', 'Sr': '锶', 'Y': '钇', 'Zr': '锆', 'Nb': '铌', 'Mo': '钼', 'Tc': '锝', 'Ru': '钌', 'Rh': '铑', 'Pd': '钯', 'Ag': '银', 'Cd': '镉',
            'In': '铟', 'Sn': '锡', 'Sb': '锑', 'Te': '碲', 'I': '碘', 'Xe': '氙'
        };

        // 清空现有选项
        select.innerHTML = '';

        // 获取并排序原子列表
        const atoms = Object.entries(window.SlaterBasis)
            .filter(([key, val]) => val.Z) // 过滤元数据
            .sort((a, b) => a[1].Z - b[1].Z);

        atoms.forEach(([symbol, data]) => {
            const option = document.createElement('option');
            option.value = symbol;
            const cnName = cnNames[symbol] || data.name;
            option.textContent = `${symbol} - ${cnName}`;
            if (symbol === 'H') option.selected = true;
            select.appendChild(option);
        });

        // 触发一次更新以确保UI一致
        // select.dispatchEvent(new Event('change')); // 不自动触发，保持默认 H
    }

    // 更新原子选择器禁用状态的函数
    window.ElectronCloud.UI.updateAtomSelectorState = function () {
        const state = window.ElectronCloud.state;
        const atomSettingsBtn = document.getElementById('atom-settings-btn');
        const atomFeatureBox = document.getElementById('atom-feature-box');
        const isDisabled = state.isDrawing || state.pointCount > 0;

        if (atomSettingsBtn) {
            atomSettingsBtn.disabled = isDisabled;
            atomSettingsBtn.style.opacity = isDisabled ? '0.5' : '1';
        }
        if (atomFeatureBox) {
            atomFeatureBox.style.opacity = isDisabled ? '0.5' : '1';
            atomFeatureBox.style.pointerEvents = isDisabled ? 'none' : 'auto';
        }
    };

    /**
     * 从轨道键中提取基础轨道类型（如 '2p_x' -> '2p', '3d_xy' -> '3d'）
     * @param {string} orbitalKey - 轨道键
     * @returns {string} - 基础轨道键
     */
    window.ElectronCloud.UI.getBaseOrbitalKey = function (orbitalKey) {
        if (!orbitalKey) return '';
        // 匹配 'nl' 部分，如 '1s', '2p', '3d', '4f' 等
        const match = orbitalKey.match(/^(\d+[spdf])/);
        return match ? match[1] : orbitalKey;
    };

    /**
     * 检测原子是否支持指定轨道
     * @param {string} atomSymbol - 原子符号（如 'H', 'C', 'Fe'）
     * @param {string} orbitalKey - 轨道键（如 '2p_x', '3d_xy'）
     * @returns {boolean} - 是否支持
     */
    window.ElectronCloud.UI.atomHasOrbital = function (atomSymbol, orbitalKey) {
        if (!window.SlaterBasis || !window.SlaterBasis[atomSymbol]) {
            // 对于H原子，特殊处理：只支持类氢轨道（任意n,l满足l<n）
            if (atomSymbol === 'H') {
                return true; // 氢原子支持所有轨道（类氢模型）
            }
            return false;
        }

        const atomData = window.SlaterBasis[atomSymbol];
        const orbitals = atomData.orbitals;
        if (!orbitals) return false;

        // 获取基础轨道键
        const baseKey = window.ElectronCloud.UI.getBaseOrbitalKey(orbitalKey);

        // 检查原子是否有该轨道
        return baseKey in orbitals;
    };

    /**
     * 更新原子下拉列表中各选项的禁用状态
     * 根据当前选中的轨道，禁用不支持这些轨道的原子
     */
    window.ElectronCloud.UI.updateAtomOptionsDisableState = function () {
        const state = window.ElectronCloud.state;
        const ui = window.ElectronCloud.ui;
        const atomSelect = document.getElementById('atom-select');

        if (!atomSelect || !window.SlaterBasis) return;

        // 获取当前选中的轨道列表
        const selectedOrbitals = state.currentOrbitals && state.currentOrbitals.length > 0
            ? state.currentOrbitals
            : [state.currentOrbital || '1s'];

        // 检查是否在比照模式（比照模式下不禁用，因为每个slot有自己的原子选择）
        const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
        if (isCompareMode) {
            // 比照模式下恢复所有原子可选
            Array.from(atomSelect.options).forEach(option => {
                option.disabled = false;
                option.style.opacity = '';
                option.style.color = '';
            });
            return;
        }

        // 遍历所有原子选项
        Array.from(atomSelect.options).forEach(option => {
            const atomSymbol = option.value;

            // 检查该原子是否支持所有选中的轨道
            let supportsAll = true;
            for (const orbitalKey of selectedOrbitals) {
                if (!window.ElectronCloud.UI.atomHasOrbital(atomSymbol, orbitalKey)) {
                    supportsAll = false;
                    break;
                }
            }

            // 更新禁用状态
            if (supportsAll) {
                option.disabled = false;
                option.style.opacity = '';
                option.style.color = '';
            } else {
                option.disabled = true;
                option.style.opacity = '0.4';
                option.style.color = '#666';
            }
        });

        // 如果当前选中的原子被禁用了，自动切换到第一个可用的原子
        const currentOption = atomSelect.options[atomSelect.selectedIndex];
        if (currentOption && currentOption.disabled) {
            // 找到第一个可用的原子
            for (const option of atomSelect.options) {
                if (!option.disabled) {
                    atomSelect.value = option.value;
                    // 触发change事件更新状态
                    atomSelect.dispatchEvent(new Event('change', { bubbles: true }));
                    break;
                }
            }
        }
    };

};

// 锁定摄像机到指定轴视角
window.ElectronCloud.UI.lockCameraToAxis = function (axis) {
    const state = window.ElectronCloud.state;
    const distance = state.camera.position.length() || 100;

    // 清除之前的自定义旋转处理
    window.ElectronCloud.UI.clearLockedAxisRotation();

    // 重置所有场景对象的旋转到初始状态（使用公共函数）
    window.ElectronCloud.UI.resetSceneObjectsRotation();

    // 重置自动旋转累计角度
    if (state.autoRotate) {
        state.autoRotate.totalAngle = 0;
    }

    // 重置所有轴的可见性，然后隐藏被锁定的轴
    if (window.ElectronCloud.Scene.resetAxesVisibility) {
        window.ElectronCloud.Scene.resetAxesVisibility();
    }
    if (window.ElectronCloud.Scene.setAxisVisibility) {
        window.ElectronCloud.Scene.setAxisVisibility(axis, false);
    }

    // 锁定哪个轴就是哪个轴指向屏幕外侧（指向观察者）
    // 先重置相机up向量为默认值
    state.camera.up.set(0, 1, 0);

    switch (axis) {
        case 'x':
            // X轴指向屏幕外侧（指向观察者），相机在+X方向看向原点
            state.camera.position.set(distance, 0, 0);
            break;
        case 'y':
            // Y轴指向屏幕外侧（指向观察者），相机在+Y方向看向原点
            state.camera.position.set(0, distance, 0);
            // Y轴向上时需要调整up向量避免万向锁
            state.camera.up.set(0, 0, -1);
            break;
        case 'z':
            // Z轴指向屏幕外侧（指向观察者），相机在+Z方向看向原点
            state.camera.position.set(0, 0, distance);
            break;
    }

    state.camera.lookAt(0, 0, 0);
    state.controls.target.set(0, 0, 0);

    // 记录当前锁定的轴
    state.lockedAxis = axis;

    // 设置自定义旋转处理（绕锁定轴旋转场景对象）
    // 注意：setupLockedAxisRotation 会禁用 OrbitControls
    window.ElectronCloud.UI.setupLockedAxisRotation(axis);
};

// 设置锁定轴的自定义旋转处理
window.ElectronCloud.UI.setupLockedAxisRotation = function (axis) {
    const state = window.ElectronCloud.state;
    const canvas = state.renderer.domElement;

    // 先清除之前的监听器
    window.ElectronCloud.UI.clearLockedAxisRotation();

    // **关键修复**：完全禁用 OrbitControls 以避免事件冲突
    // OrbitControls 即使 enableRotate=false，仍会拦截 pointer 事件
    state.controls.enabled = false;

    let isDragging = false;
    let previousMouseX = 0;
    let previousMouseY = 0;

    // 缩放相关
    let initialPinchDistance = 0;
    let isPinching = false;

    // 获取旋转轴向量
    const getRotationAxis = function () {
        switch (state.lockedAxis) {
            case 'x': return new THREE.Vector3(1, 0, 0);
            case 'y': return new THREE.Vector3(0, 1, 0);
            case 'z': return new THREE.Vector3(0, 0, 1);
            default: return null;
        }
    };

    // 旋转场景对象
    const rotateSceneObjects = function (angle) {
        const rotationAxis = getRotationAxis();
        if (!rotationAxis) return;

        if (state.points) {
            state.points.rotateOnWorldAxis(rotationAxis, angle);
        }
        if (state.angularOverlay) {
            state.angularOverlay.rotateOnWorldAxis(rotationAxis, angle);
        }
        if (state.customAxes) {
            state.customAxes.rotateOnWorldAxis(rotationAxis, angle);
        }
    };

    // 缩放处理
    const handleZoom = function (delta) {
        const zoomSpeed = 0.001;
        const factor = 1 - delta * zoomSpeed;
        const minDistance = 5;
        const maxDistance = 500;

        const currentDistance = state.camera.position.length();
        const newDistance = Math.max(minDistance, Math.min(maxDistance, currentDistance * factor));

        state.camera.position.normalize().multiplyScalar(newDistance);
    };

    // 使用 pointer 事件以获得更好的兼容性（包括触摸）
    const onPointerDown = function (e) {
        if (!state.lockedAxis) return;

        // 左键或触摸开始旋转
        if (e.button === 0 || e.pointerType === 'touch') {
            isDragging = true;
            previousMouseX = e.clientX;
            previousMouseY = e.clientY;
            canvas.setPointerCapture(e.pointerId);
        }

        e.preventDefault();
    };

    const onPointerMove = function (e) {
        if (!state.lockedAxis) return;
        if (!isDragging) return;

        const deltaX = e.clientX - previousMouseX;
        // deltaY 暂时不用，但保留以备将来使用
        // const deltaY = e.clientY - previousMouseY;

        previousMouseX = e.clientX;
        previousMouseY = e.clientY;

        // 旋转速度
        const rotationSpeed = 0.01;
        const angle = deltaX * rotationSpeed;

        rotateSceneObjects(angle);

        e.preventDefault();
    };

    const onPointerUp = function (e) {
        if (isDragging) {
            isDragging = false;
            try {
                canvas.releasePointerCapture(e.pointerId);
            } catch (err) {
                // 忽略释放捕获的错误
            }
        }
    };

    // 滚轮缩放（包括触控板双指缩放）
    const onWheel = function (e) {
        if (!state.lockedAxis) return;

        // 触控板双指缩放通常会设置 ctrlKey，并且 deltaY 值较小
        // 普通滚轮的 deltaY 通常是 100 或 -100 的倍数
        if (e.ctrlKey) {
            // 触控板捏合缩放：deltaY > 0 表示捏合（放大手势），应该缩小视图
            // 但用户习惯是捏合=缩小，所以需要反转
            handleZoom(-e.deltaY * 5);
        } else {
            // 普通滚轮
            handleZoom(e.deltaY);
        }
        e.preventDefault();
    };

    // 触摸缩放（双指捏合）
    const activeTouches = new Map();

    const onTouchStart = function (e) {
        for (let touch of e.changedTouches) {
            activeTouches.set(touch.identifier, { x: touch.clientX, y: touch.clientY });
        }

        if (activeTouches.size === 2) {
            const touches = Array.from(activeTouches.values());
            initialPinchDistance = Math.hypot(
                touches[0].x - touches[1].x,
                touches[0].y - touches[1].y
            );
            isPinching = true;
        }
    };

    const onTouchMove = function (e) {
        for (let touch of e.changedTouches) {
            if (activeTouches.has(touch.identifier)) {
                activeTouches.set(touch.identifier, { x: touch.clientX, y: touch.clientY });
            }
        }

        if (isPinching && activeTouches.size === 2) {
            const touches = Array.from(activeTouches.values());
            const currentDistance = Math.hypot(
                touches[0].x - touches[1].x,
                touches[0].y - touches[1].y
            );

            const delta = (initialPinchDistance - currentDistance) * 2;
            handleZoom(delta);
            initialPinchDistance = currentDistance;

            e.preventDefault();
        }
    };

    const onTouchEnd = function (e) {
        for (let touch of e.changedTouches) {
            activeTouches.delete(touch.identifier);
        }

        if (activeTouches.size < 2) {
            isPinching = false;
        }
    };

    // 保存引用以便后续移除
    state.lockedAxisHandlers = {
        canvas: canvas,
        pointerdown: onPointerDown,
        pointermove: onPointerMove,
        pointerup: onPointerUp,
        wheel: onWheel,
        touchstart: onTouchStart,
        touchmove: onTouchMove,
        touchend: onTouchEnd
    };

    // 绑定事件
    canvas.addEventListener('pointerdown', onPointerDown);
    canvas.addEventListener('pointermove', onPointerMove);
    canvas.addEventListener('pointerup', onPointerUp);
    canvas.addEventListener('pointercancel', onPointerUp);
    canvas.addEventListener('wheel', onWheel, { passive: false });
    canvas.addEventListener('touchstart', onTouchStart, { passive: true });
    canvas.addEventListener('touchmove', onTouchMove, { passive: false });
    canvas.addEventListener('touchend', onTouchEnd, { passive: true });

    // 也在 window 上监听 pointerup 以确保拖动到窗口外也能释放
    window.addEventListener('pointerup', onPointerUp);
};

// 清除锁定轴的自定义旋转处理
window.ElectronCloud.UI.clearLockedAxisRotation = function () {
    const state = window.ElectronCloud.state;

    if (state.lockedAxisHandlers) {
        const handlers = state.lockedAxisHandlers;
        const canvas = handlers.canvas;

        if (canvas) {
            // 移除 pointer 事件
            if (handlers.pointerdown) {
                canvas.removeEventListener('pointerdown', handlers.pointerdown);
                canvas.removeEventListener('pointermove', handlers.pointermove);
                canvas.removeEventListener('pointerup', handlers.pointerup);
                canvas.removeEventListener('pointercancel', handlers.pointerup);
            }
            // 移除滚轮事件
            if (handlers.wheel) {
                canvas.removeEventListener('wheel', handlers.wheel);
            }
            // 移除触摸事件
            if (handlers.touchstart) {
                canvas.removeEventListener('touchstart', handlers.touchstart);
                canvas.removeEventListener('touchmove', handlers.touchmove);
                canvas.removeEventListener('touchend', handlers.touchend);
            }
            // 兼容旧版 mouse 事件（如果有）
            if (handlers.mousedown) {
                canvas.removeEventListener('mousedown', handlers.mousedown, true);
                canvas.removeEventListener('mousemove', handlers.mousemove, true);
                canvas.removeEventListener('mouseup', handlers.mouseup, true);
            }
        }

        // 移除 window 上的监听器
        if (handlers.pointerup) {
            window.removeEventListener('pointerup', handlers.pointerup);
        }
        if (handlers.mouseup) {
            window.removeEventListener('mouseup', handlers.mouseup, true);
        }

        state.lockedAxisHandlers = null;
    }

    // 重新启用 OrbitControls
    if (state.controls) {
        state.controls.enabled = true;
    }
};

// 自动旋转与锁定视角互斥逻辑
window.ElectronCloud.UI.updateRotationLockMutualExclusion = function () {
    const state = window.ElectronCloud.state;
    const viewLockBtns = document.querySelectorAll('.view-lock-btn');
    const rotationSubmenu = document.getElementById('rotation-submenu');
    const rotationSettingsBtn = document.getElementById('rotation-settings-btn');
    const rotationSpeedRange = document.getElementById('rotation-speed-range');
    const rotationAxisInputs = document.querySelectorAll('.axis-inputs input');
    const autoRotateToggle = document.getElementById('auto-rotate-toggle');
    const rotationFeatureBox = document.getElementById('rotation-feature-box');

    const isAutoRotating = state.autoRotate && state.autoRotate.enabled;
    const isViewLocked = state.lockedAxis !== null && state.lockedAxis !== undefined;

    // 如果自动旋转启用，禁用锁定视角按钮
    viewLockBtns.forEach(btn => {
        if (isAutoRotating) {
            btn.disabled = true;
            btn.style.opacity = '0.4';
            btn.style.pointerEvents = 'none';
        } else {
            btn.disabled = false;
            btn.style.opacity = '';
            btn.style.pointerEvents = '';
        }
    });

    // 如果视角锁定，禁用自动旋转控件
    if (autoRotateToggle) {
        if (isViewLocked) {
            autoRotateToggle.disabled = true;
            if (rotationFeatureBox) {
                rotationFeatureBox.style.opacity = '0.4';
                rotationFeatureBox.style.pointerEvents = 'none';
            }
        } else {
            autoRotateToggle.disabled = false;
            if (rotationFeatureBox) {
                rotationFeatureBox.style.opacity = '';
                rotationFeatureBox.style.pointerEvents = '';
            }
        }
    }

    if (rotationSpeedRange) {
        if (isViewLocked) {
            rotationSpeedRange.disabled = true;
            rotationSpeedRange.style.opacity = '0.4';
        } else {
            rotationSpeedRange.disabled = false;
            rotationSpeedRange.style.opacity = '';
        }
    }

    rotationAxisInputs.forEach(input => {
        if (isViewLocked) {
            input.disabled = true;
            input.style.opacity = '0.4';
        } else {
            input.disabled = false;
            input.style.opacity = '';
        }
    });

    if (rotationSettingsBtn) {
        if (isViewLocked) {
            rotationSettingsBtn.disabled = true;
            rotationSettingsBtn.style.opacity = '0.4';
        } else {
            rotationSettingsBtn.disabled = false;
            rotationSettingsBtn.style.opacity = '';
        }
    }
};

// ========================================
// 公共工具函数（避免代码重复）
// ========================================

// 计算坐标轴缩放比例（公共函数，避免重复代码）
window.ElectronCloud.UI.calculateAxesScale = function (scaleFactor) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;

    if (scaleFactor === 0 || state.farthestDistance === 0) {
        return { visible: false, scale: 1 };
    }

    const orbitalRadius = Math.max(constants.AXES_BASE_SIZE, state.farthestDistance);
    const targetSize = orbitalRadius * scaleFactor;
    const scale = targetSize / constants.AXES_BASE_SIZE;

    return { visible: true, scale: scale };
};

// 清除轨道选择框的所有选中状态（公共函数，避免重复代码）
window.ElectronCloud.UI.resetOrbitalSelections = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    if (!ui.orbitalSelect) return;

    // 清除所有选中状态
    Array.from(ui.orbitalSelect.options).forEach(option => {
        option.selected = false;
        option.classList.remove('force-selected', 'force-unselected', 'compare-color-0', 'compare-color-1', 'compare-color-2');
    });

    state.currentOrbitals = [];

    // 同步更新自定义列表视觉状态
    window.ElectronCloud.UI.updateCustomListVisuals();

    // 更新选择计数
    window.ElectronCloud.UI.updateSelectionCount();

    // 【新增】重置比照模式选择器
    if (ui.compareToggle && ui.compareToggle.checked) {
        window.ElectronCloud.UI.resetCompareOrbitalSelectors();
    }
};

// 重置所有场景对象的旋转状态（公共函数，避免重复代码）
// 仅重置旋转，不清除锁定状态或其他逻辑
window.ElectronCloud.UI.resetSceneObjectsRotation = function () {
    const state = window.ElectronCloud.state;

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
};

// ========================================

// 处理坐标系大小调节
window.ElectronCloud.UI.onAxesSizeChange = function (event) {
    const state = window.ElectronCloud.state;

    if (state.customAxes) {
        const value = parseInt(event.target.value, 10);

        // 将滑动条值转换为比例系数 (0-1.05)
        const scaleFactor = value / 100;
        state.axesScaleFactor = scaleFactor;

        // 使用公共函数计算缩放
        const result = window.ElectronCloud.UI.calculateAxesScale(scaleFactor);
        state.customAxes.visible = result.visible;

        if (result.visible) {
            state.customAxes.scale.set(result.scale, result.scale, result.scale);
        }
    }
};

// 更新坐标轴滑动条的可用状态
window.ElectronCloud.UI.updateAxesSizeRangeState = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    if (!ui.axesSizeRange) return;

    // 渲染中时置灰不可修改，其他时候可修改
    const isDisabled = state.isDrawing && !state.samplingCompleted;
    ui.axesSizeRange.disabled = isDisabled;

    const label = document.querySelector('label[for="axes-size-range"]');
    if (label) {
        label.style.opacity = isDisabled ? '0.5' : '1';
        label.title = isDisabled ? '渲染中不可修改' : '';
    }
};

// 处理角向分布叠加显隐
window.ElectronCloud.UI.onAngularOverlayToggle = function (event) {
    const state = window.ElectronCloud.state;

    if (event.target.checked) {
        // 检查是否有采样数据或采样是否完成
        if (state.pointCount === 0) {
            event.target.checked = false;
            alert('请先开始采样后再启用角向分布3D显示');
            return;
        }

        if (state.isDrawing && state.pointCount < state.MAX_POINTS) {
            event.target.checked = false;
            alert('请等待采样完成后再启用角向分布3D显示');
            return;
        }

        window.ElectronCloud.Visualization.updateAngularOverlay();
    } else {
        window.ElectronCloud.Scene.clearAngularOverlay();
    }
};

// 处理95%等值面轮廓显隐（网格方案）
window.ElectronCloud.UI.onContourOverlayToggle = function (event) {
    const state = window.ElectronCloud.state;

    if (event.target.checked) {
        // 检查是否有采样数据
        if (state.pointCount === 0) {
            event.target.checked = false;
            alert('请先开始采样后再启用轨道轮廓显示');
            return;
        }

        if (state.isDrawing && state.pointCount < state.MAX_POINTS) {
            event.target.checked = false;
            alert('请等待采样完成后再启用轨道轮廓显示');
            return;
        }

        // 使用网格等值面方案
        window.ElectronCloud.Visualization.updateContourOverlay();
    } else {
        // 移除等值面网格
        if (state.contourOverlay) {
            state.scene.remove(state.contourOverlay);
            state.contourOverlay.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
            state.contourOverlay = null;
        }
    }
};

// 处理轨道选择变化
window.ElectronCloud.UI.onOrbitalSelectChange = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    const opts = Array.from(ui.orbitalSelect.selectedOptions || []).map(o => o.value);

    // 所有模式下统一处理：有选中则使用选中的，没有则为空数组
    // 【修改】单选模式现在也不再默认选择 '1s'，与多选/比照模式保持一致
    state.currentOrbitals = opts.length ? opts : [];

    // 单选模式下手动更新选中样式（因为浏览器不会自动更新 option:checked 背景）
    if (!(ui.multiselectToggle && ui.multiselectToggle.checked) && !(ui.compareToggle && ui.compareToggle.checked)) {
        window.ElectronCloud.UI.updateSingleSelectStyle();
    }

    // 重新选择轨道后，重置3D角向开关为关闭状态
    if (ui.angular3dToggle && ui.angular3dToggle.checked) {
        ui.angular3dToggle.checked = false;
        const changeEvent = new Event('change', { bubbles: true });
        ui.angular3dToggle.dispatchEvent(changeEvent);
    }

    // 重新选择轨道后，重置轨道轮廓开关为关闭状态
    if (ui.contour3dToggle && ui.contour3dToggle.checked) {
        ui.contour3dToggle.checked = false;
        const changeEvent = new Event('change', { bubbles: true });
        ui.contour3dToggle.dispatchEvent(changeEvent);
    }

    // 更新多选模式或比照模式的选择计数
    if ((ui.multiselectToggle && ui.multiselectToggle.checked) || (ui.compareToggle && ui.compareToggle.checked)) {
        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();

        // 更新清空按钮状态
        window.ElectronCloud.UI.updateClearAllSelectionsState();
    }

    // 【新增】更新原子选项禁用状态（根据选中轨道禁用不支持的原子）
    if (window.ElectronCloud.UI.updateAtomOptionsDisableState) {
        window.ElectronCloud.UI.updateAtomOptionsDisableState();
    }
};

// 单选模式下更新选中样式
window.ElectronCloud.UI.updateSingleSelectStyle = function () {
    const ui = window.ElectronCloud.ui;
    if (!ui.orbitalSelect) return;

    // 清除所有 option 的高亮样式
    Array.from(ui.orbitalSelect.options).forEach(opt => {
        opt.classList.remove('selected-green');
        opt.style.background = '';
        opt.style.color = '';
        opt.style.fontWeight = '';
    });

    // 给当前选中的 option 添加高亮样式
    if (ui.orbitalSelect.selectedIndex >= 0) {
        const selectedOpt = ui.orbitalSelect.options[ui.orbitalSelect.selectedIndex];
        selectedOpt.classList.add('selected-green');
        // 强制内联样式，确保在所有浏览器中生效
        selectedOpt.style.background = 'rgba(80, 200, 120, 0.35)';
        selectedOpt.style.color = 'rgba(120, 220, 150, 1)';
        selectedOpt.style.fontWeight = 'bold';
    }
};

// 更新角向分布3D开关的可用状态
window.ElectronCloud.UI.updateAngular3DToggleState = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!ui.angular3dToggle) return;

    const angular3dLabel = document.querySelector('label[for="angular-3d-toggle"]');

    // 只在对比模式下直接禁用，不管其他条件
    if (ui.compareToggle && ui.compareToggle.checked) {
        ui.angular3dToggle.disabled = true;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = '对比模式下不可用';
        }
        // 【关键修复】比照模式下虽然禁用了角向形状，但必须更新（启用）轨道轮廓开关
        window.ElectronCloud.UI.updateContour3DToggleState();
        return;
    }

    // 【新增】多选模式禁用3D角向形状
    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        ui.angular3dToggle.disabled = true;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = '多选模式下不可用';
        }
        // 同步更新轨道轮廓状态
        window.ElectronCloud.UI.updateContour3DToggleState();
        return;
    }

    // 【杂化模式】禁用3D角向形状
    if (state.isHybridMode === true) {
        ui.angular3dToggle.disabled = true;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = '杂化模式下不可用';
        }
        // 但仍需更新轨道轮廓状态（轨道轮廓在杂化模式下应可用）
        window.ElectronCloud.UI.updateContour3DToggleState();
        return;
    }

    // 在采样过程中禁用开关
    if (state.isDrawing && state.pointCount < state.MAX_POINTS) {
        ui.angular3dToggle.disabled = true;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = '采样完成后可用';
        }
    } else if (state.pointCount === 0) {
        ui.angular3dToggle.disabled = true;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = '请先开始采样';
        }
    } else {
        ui.angular3dToggle.disabled = false;
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '1';
            angular3dLabel.style.cursor = 'pointer';
            angular3dLabel.title = '';
        }
    }

    // 同步更新轨道轮廓开关（保持一致的逻辑）
    window.ElectronCloud.UI.updateContour3DToggleState();
};

// 更新轨道轮廓开关的可用状态（与3D角向形状保持一致）
window.ElectronCloud.UI.updateContour3DToggleState = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    const contour3dToggle = document.getElementById('contour-3d-toggle');
    if (!contour3dToggle) return;

    const contour3dLabel = document.querySelector('label[for="contour-3d-toggle"]');

    // 【修改】比照模式下不禁用轨道轮廓，允许多色显示

    // 在采样过程中禁用
    if (state.isDrawing && state.pointCount < state.MAX_POINTS) {
        contour3dToggle.disabled = true;
        if (contour3dLabel) {
            contour3dLabel.style.opacity = '0.5';
            contour3dLabel.style.cursor = 'not-allowed';
            contour3dLabel.title = '采样完成后可用';
        }
    } else if (state.pointCount === 0) {
        contour3dToggle.disabled = true;
        if (contour3dLabel) {
            contour3dLabel.style.opacity = '0.5';
            contour3dLabel.style.cursor = 'not-allowed';
            contour3dLabel.title = '请先开始采样';
        }
    } else {
        contour3dToggle.disabled = false;
        if (contour3dLabel) {
            contour3dLabel.style.opacity = '1';
            contour3dLabel.style.cursor = 'pointer';
            contour3dLabel.title = '';
        }
    }
};

// 控件辅助函数
window.ElectronCloud.UI.onSpeedChange = function () {
    // 速度变化直接由采样模块读取
};

// 轨道预估最远距离的先验知识表（基于量子数n）
// 公式约为 15 * n^2，这里直接查表避免实时计算
window.ElectronCloud.UI.orbitalDistanceTable = {
    // n=1 轨道，预估最远距离约 15
    '1s': 15,
    // n=2 轨道，预估最远距离约 60
    '2s': 60, '2px': 60, '2py': 60, '2pz': 60,
    // n=3 轨道，预估最远距离约 135
    '3s': 135, '3px': 135, '3py': 135, '3pz': 135,
    '3d_xy': 135, '3d_xz': 135, '3d_yz': 135, '3d_x2-y2': 135, '3d_z2': 135,
    // n=4 轨道，预估最远距离约 240
    '4s': 240, '4px': 240, '4py': 240, '4pz': 240,
    '4d_xy': 240, '4d_xz': 240, '4d_yz': 240, '4d_x2-y2': 240, '4d_z2': 240,
    '4f_xyz': 240, '4f_xz2': 240, '4f_yz2': 240, '4f_z3': 240,
    '4f_z(x2-y2)': 240, '4f_x(x2-3y2)': 240, '4f_y(3x2-y2)': 240
};

// 获取点大小（根据轨道尺寸等比例缩放）
// 以"最远距离 60"为基准（对应 n=2 轨道）
window.ElectronCloud.UI.getPointSize = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;
    const Hydrogen = window.Hydrogen;
    const distanceTable = window.ElectronCloud.UI.orbitalDistanceTable;

    if (!ui.sizeRange) return 0.06;

    // 滑条范围 1-200，默认值 50
    // 原来的中等大小 0.06 对应滑条值 50
    // 滑条值 100 对应原来的大小 0.12，滑条值 200 对应 0.24
    const sliderValue = parseInt(ui.sizeRange.value);
    // 线性映射：1->0.012, 50->0.06, 100->0.12, 200->0.24
    const baseSize = (sliderValue / 50) * 0.06;

    // 基准最远距离（对应 n=2 轨道）
    const baseDistance = 60;
    const atomType = state.currentAtom || 'H';

    // 获取当前选中轨道中最大的预估距离
    let maxDistance = baseDistance; // 默认值

    // 优先使用动态估算函数（支持非氢原子）
    if (Hydrogen && Hydrogen.estimateOrbitalRadius95 && state.currentOrbitals && state.currentOrbitals.length > 0) {
        for (const orbitalKey of state.currentOrbitals) {
            const dist = Hydrogen.estimateOrbitalRadius95(atomType, orbitalKey);
            if (dist && dist > maxDistance) {
                maxDistance = dist;
            }
        }
    }
    // 备用：从静态表获取（仅氢原子）
    else if (state.currentOrbitals && state.currentOrbitals.length > 0) {
        for (const orbitalKey of state.currentOrbitals) {
            const dist = distanceTable[orbitalKey];
            if (dist && dist > maxDistance) {
                maxDistance = dist;
            }
        }
    }
    // 再备用：从 UI 选择框获取
    else if (ui.orbitalSelect) {
        const selectedOptions = ui.orbitalSelect.selectedOptions;
        if (selectedOptions && selectedOptions.length > 0) {
            for (const option of selectedOptions) {
                const dist = distanceTable[option.value];
                if (dist && dist > maxDistance) {
                    maxDistance = dist;
                }
            }
        } else if (ui.orbitalSelect.value) {
            // 单选模式
            const dist = distanceTable[ui.orbitalSelect.value];
            if (dist) {
                maxDistance = dist;
            }
        }
    }

    // 点大小与轨道尺寸等比例缩放
    const scaleFactor = maxDistance / baseDistance;
    return baseSize * scaleFactor;
};

window.ElectronCloud.UI.onSizeChange = function () {
    const state = window.ElectronCloud.state;

    if (state.points && state.points.material) {
        state.points.material.size = window.ElectronCloud.UI.getPointSize();
    }
};

window.ElectronCloud.UI.onOpacityChange = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!ui.opacityRange) return;

    const value = parseInt(ui.opacityRange.value, 10); // 0-100

    // 新映射：前50%映射到透明度 (0-50 -> 0.05-1.0)
    // 后50%映射到辉光强度 (50-100 -> 0-4.5)，使用曲线映射
    // 默认值为67（约2/3位置），处于辉光区间

    if (value <= 50) {
        // 透明度区间：0-50 映射到 0.05-1.0
        const opacity = 0.05 + (value / 50) * 0.95;

        if (state.points && state.points.material) {
            state.points.material.opacity = opacity;
        }

        // 关闭辉光
        if (state.bloomPass) {
            state.bloomPass.strength = 0;
        }
        state.bloomEnabled = false;

    } else {
        // 辉光区间：50-100 映射到辉光强度 0-4.5
        // 使用曲线映射(指数函数)，在低值时能精细调节
        // 透明度保持最大
        if (state.points && state.points.material) {
            state.points.material.opacity = 1.0;
        }

        // 开启辉光，使用曲线映射计算强度
        const glowProgress = (value - 50) / 50; // 0-1
        // 使用平方曲线: y = x^2，这样小值区间更精细
        const curvedProgress = glowProgress * glowProgress;
        const bloomStrength = curvedProgress * 4.5; // 0-4.5

        if (state.bloomPass) {
            state.bloomPass.strength = bloomStrength;
        }
        state.bloomEnabled = true;
        state.bloomStrength = bloomStrength;
    }
};

window.ElectronCloud.UI.onCenterLockChange = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!state.controls || !ui.centerLock) return;

    const locked = ui.centerLock.checked;
    // 更新 center-lock-box 的 active 样式
    const centerLockBox = document.getElementById('center-lock-box');
    if (centerLockBox) {
        centerLockBox.classList.toggle('active', locked);
    }

    // 禁用平移并锁定目标在原点
    state.controls.enablePan = !locked;
    if (locked) {
        state.controls.target.set(0, 0, 0);
        state.controls.update();
    }
};

// 处理相位显示开关
window.ElectronCloud.UI.onPhaseToggle = function (event) {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    state.usePhaseColoring = event.target.checked;

    // 更新 phase-box 的 active 样式
    const phaseBox = document.getElementById('phase-box');
    if (phaseBox) {
        phaseBox.classList.toggle('active', state.usePhaseColoring);
    }

    // 只有在有点云数据时才更新颜色
    if (state.points && state.pointCount > 0) {
        // 更新点云颜色
        if (window.ElectronCloud.Visualization.updatePointColors) {
            window.ElectronCloud.Visualization.updatePointColors();
        }

        // 【关键修复】同步更新 baseColors 以确保闪烁模式使用正确的相位颜色
        if (state.baseColors && state.points.geometry) {
            const colors = state.points.geometry.attributes.color;
            if (colors) {
                const pointCount = state.pointCount || 0;
                for (let i = 0; i < pointCount * 3; i++) {
                    state.baseColors[i] = colors.array[i];
                }
                state.baseColorsCount = pointCount;
            }
        }

        // 更新轮廓颜色（如果轮廓已显示）
        if (state.contourOverlay && window.ElectronCloud.Visualization.updateContourOverlay) {
            // 重新生成等值面以应用新颜色
            window.ElectronCloud.Visualization.updateContourOverlay();
        } else if (state.hybridContourOverlay && window.ElectronCloud.Visualization.createHybridContourOverlays) {
            // 重新生成杂化轨道等值面
            window.ElectronCloud.Visualization.createHybridContourOverlays();
        }
    }
};

// 辅助函数：检查并同步 centerLock 状态
// 用于其他模块恢复状态时调用
window.ElectronCloud.UI.syncCenterLockState = function () {
    const state = window.ElectronCloud.state;
    if (!state.controls) return;

    const centerLock = document.getElementById('center-lock');
    const locked = centerLock && centerLock.checked;

    state.controls.enablePan = !locked;
    if (locked) {
        state.controls.target.set(0, 0, 0);
    }
};

// 渲染相位切换
window.ElectronCloud.UI.onPhaseToggleChange = function () {
    const phaseToggle = document.getElementById('phase-toggle');
    const phaseBox = document.getElementById('phase-box');

    if (phaseBox && phaseToggle) {
        phaseBox.classList.toggle('active', phaseToggle.checked);
    }
};

// 为所有 mode-toggle-box 添加点击切换逻辑
window.ElectronCloud.UI.setupModeToggleBoxes = function () {
    const boxes = document.querySelectorAll('.mode-toggle-box');
    boxes.forEach((box, index) => {
        box.addEventListener('click', function (e) {
            // 如果点击的是 checkbox 本身，不需要额外处理
            if (e.target.tagName === 'INPUT') {
                return;
            }

            const checkbox = box.querySelector('input[type="checkbox"]');
            if (checkbox && !checkbox.disabled) {
                checkbox.checked = !checkbox.checked;
                // 触发 change 事件
                checkbox.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });

        // 设置鼠标样式为指针
        box.style.cursor = 'pointer';
    });

    // 为 phase-toggle 添加 change 监听
    const phaseToggle = document.getElementById('phase-toggle');
    if (phaseToggle) {
        phaseToggle.addEventListener('change', window.ElectronCloud.UI.onPhaseToggleChange);
    }

    // 初始化 phase-box 和 center-lock-box 的 active 状态
    window.ElectronCloud.UI.onPhaseToggleChange();
    window.ElectronCloud.UI.onCenterLockChange();
};

window.ElectronCloud.UI.onPauseToggle = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;

    if (!state.isDrawing) {
        // 继续采样，但不启用实时图表更新
        state.isDrawing = true;
        if (ui.pauseButton) ui.pauseButton.textContent = '暂停';
        window.ElectronCloud.UI.updateAngular3DToggleState();
        window.ElectronCloud.UI.updateAxesSizeRangeState();
        window.ElectronCloud.UI.updateModeSwitcherState();
        return;
    }
    // 暂停
    state.isDrawing = false;
    if (ui.pauseButton) ui.pauseButton.textContent = '继续';
    window.ElectronCloud.UI.updateAngular3DToggleState();
    window.ElectronCloud.UI.updateAxesSizeRangeState();
    window.ElectronCloud.UI.updateModeSwitcherState();
};

// 多选模式功能
window.ElectronCloud.UI.onMultiselectToggle = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    // 渲染开始后禁止切换模式
    const hasStarted = state.pointCount > 0 || state.isDrawing;
    if (hasStarted) {
        // 恢复checkbox状态
        if (ui.multiselectToggle) {
            ui.multiselectToggle.checked = !ui.multiselectToggle.checked;
        }
        return;
    }

    // 【强制清除】切换模式时，无论如何都先清除所有选择
    window.ElectronCloud.UI.resetOrbitalSelections();

    // 互斥：取消比照模式
    if (ui.compareToggle && ui.compareToggle.checked) {
        ui.compareToggle.checked = false;
        // 不递归调用onCompareToggle，避免重复清除
        const compareBox = document.getElementById('compare-box');
        if (compareBox) compareBox.classList.remove('active');

        // 【修复】从比照模式切换过来时，恢复相位开关的可用状态
        const phaseToggle = document.getElementById('phase-toggle');
        const phaseBox = document.getElementById('phase-box');
        if (phaseToggle) {
            phaseToggle.disabled = false;
        }
        if (phaseBox) {
            phaseBox.classList.remove('disabled');
        }

        // 重新启用显示开关功能
        window.ElectronCloud.UI.enableAngular3DToggle();
    }

    const isMultiselect = ui.multiselectToggle && ui.multiselectToggle.checked;
    const label = document.querySelector('label[for="orbital-select"]');
    const controlPanel = document.getElementById('control-panel');
    const multiselectControls = document.querySelector('.multiselect-controls');
    const multiselectBox = document.getElementById('multiselect-box');

    // 更新多选框激活样式
    if (multiselectBox) {
        multiselectBox.classList.toggle('active', isMultiselect);
    }

    if (isMultiselect) {
        // 【修复】明确恢复相位开关的可用状态（多选模式支持相位显示）
        const phaseToggle = document.getElementById('phase-toggle');
        const phaseBox = document.getElementById('phase-box');
        if (phaseToggle) {
            phaseToggle.disabled = false;
        }
        if (phaseBox) {
            phaseBox.classList.remove('disabled');
        }

        if (controlPanel) controlPanel.classList.add('multiselect-active');
        if (multiselectControls) multiselectControls.classList.add('visible');
        if (label) label.style.display = 'none';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = true;
            ui.orbitalSelect.size = 5;
        }

        // 多选模式下禁用角向分布曲线选项
        window.ElectronCloud.UI.disableAngularPlotOption();

        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();
        window.ElectronCloud.UI.updateCustomListVisuals();
    } else {
        // 重新启用显示开关功能（3D角向形状在多选模式下也是可用的，所以这里不需要特别处理）
        // window.ElectronCloud.UI.enableAngular3DToggle(); // 多选模式下没有禁用，所以退出时也不需要启用

        // 【新增】重新启用角向分布曲线选项
        window.ElectronCloud.UI.enableAngularPlotOption();

        if (label) {
            label.textContent = '选择轨道';
            label.style.display = '';
        }
        if (controlPanel) controlPanel.classList.remove('multiselect-active');
        if (multiselectControls) multiselectControls.classList.remove('visible');
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = false;
            ui.orbitalSelect.size = 3; // 压缩后的高度
        }

        // 清除强制样式类（包括force-unselected）
        if (ui.orbitalSelect) {
            const options = ui.orbitalSelect.querySelectorAll('option');
            options.forEach(option => {
                option.classList.remove('force-selected', 'force-unselected');
            });
        }

        // 【修复】退出多选模式时，处理选中状态
        if (ui.orbitalSelect) {
            const selectedCount = ui.orbitalSelect.selectedOptions.length;
            if (selectedCount > 1) {
                // 多个选中时，保留第一个选中的选项
                const firstSelected = ui.orbitalSelect.selectedOptions[0];
                Array.from(ui.orbitalSelect.options).forEach(option => {
                    option.selected = (option === firstSelected);
                });
            } else if (selectedCount === 0) {
                // 【修改】没有选中时保持空选择状态，与启动时行为一致
                // 不再默认选中第一个选项
            }
            // 同步state.currentOrbitals
            const selected = Array.from(ui.orbitalSelect.selectedOptions).map(o => o.value);
            state.currentOrbitals = selected;

            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);
        }
    }
    window.ElectronCloud.UI.updateAngular3DToggleState();
};

// 比照模式功能
window.ElectronCloud.UI.onCompareToggle = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    // 渲染开始后禁止切换模式
    const hasStarted = state.pointCount > 0 || state.isDrawing;
    if (hasStarted) {
        // 恢复checkbox状态
        if (ui.compareToggle) {
            ui.compareToggle.checked = !ui.compareToggle.checked;
        }
        return;
    }

    // 【强制清除】切换模式时，无论如何都先清除所有选择
    window.ElectronCloud.UI.resetOrbitalSelections();

    // 互斥：取消多选模式
    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        ui.multiselectToggle.checked = false;
        // 不递归调用onMultiselectToggle，避免重复清除
        const multiselectBox = document.getElementById('multiselect-box');
        if (multiselectBox) multiselectBox.classList.remove('active');
        const controlPanel = document.getElementById('control-panel');
        if (controlPanel) controlPanel.classList.remove('multiselect-active');
    }

    const isCompare = ui.compareToggle && ui.compareToggle.checked;
    const label = document.querySelector('label[for="orbital-select"]');
    const controlPanel = document.getElementById('control-panel');
    const multiselectControls = document.querySelector('.multiselect-controls');
    const phaseToggle = document.getElementById('phase-toggle');
    const compareBox = document.getElementById('compare-box');

    // 更新比照框激活样式
    if (compareBox) {
        compareBox.classList.toggle('active', isCompare);
    }

    if (isCompare) {
        // 【性能优化】懒加载：首次进入比照模式时才初始化选择器
        if (!window.ElectronCloud.UI._compareSelectorsInitialized) {
            console.log('[性能优化] 懒加载初始化比照模式选择器...');
            window.ElectronCloud.UI.initCompareOrbitalSelectors();
            window.ElectronCloud.UI._compareSelectorsInitialized = true;
        }

        // 禁用渲染相位
        if (phaseToggle) {
            phaseToggle.disabled = true;
            phaseToggle.checked = false;
        }
        // 更新 phase-box 的 active 样式和禁用样式
        const phaseBox = document.getElementById('phase-box');
        if (phaseBox) {
            phaseBox.classList.add('disabled');
            phaseBox.classList.remove('active');
        }

        // 禁用3D角向形状显示，避免造成bug
        window.ElectronCloud.UI.disableAngular3DToggle();

        // 【新增】比照模式下启用角向分布曲线选项（与径向分布同步支持）
        window.ElectronCloud.UI.enableAngularPlotOption();

        // 【新增】比照模式下禁用实验面板的原子选择
        const atomSubmenu = document.getElementById('atom-submenu');
        const atomSelectMain = document.getElementById('atom-select');
        if (atomSubmenu) {
            atomSubmenu.classList.add('disabled');
            atomSubmenu.style.pointerEvents = 'none';
            atomSubmenu.style.opacity = '0.5';
            atomSubmenu.title = '比照模式下请在左侧选择器中单独设置每个轨道的原子';
        }
        if (atomSelectMain) {
            atomSelectMain.disabled = true;
        }

        if (controlPanel) {
            controlPanel.classList.add('multiselect-active');
            controlPanel.classList.add('compare-active');
        }
        // 隐藏多选标签，显示新的比照选择器
        if (multiselectControls) multiselectControls.classList.remove('visible');
        if (label) label.style.display = 'none';

        // 【新增】显示比照模式专用选择器，隐藏原有列表
        window.ElectronCloud.UI.toggleCompareOrbitalSelectors(true);

        // 重置比照选择器状态
        window.ElectronCloud.UI.resetCompareOrbitalSelectors();

        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = true;
            ui.orbitalSelect.size = 5;
        }

        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();

        // 更新清空按钮状态
        window.ElectronCloud.UI.updateClearAllSelectionsState();
    } else {
        // 重新启用渲染相位
        if (phaseToggle) {
            phaseToggle.disabled = false;
        }
        // 更新 phase-box 的样式（恢复启用状态的样式）
        const phaseBox = document.getElementById('phase-box');
        if (phaseBox) {
            phaseBox.classList.remove('disabled');
        }

        // 重新启用3D角向形状显示
        window.ElectronCloud.UI.enableAngular3DToggle();

        // 【新增】退出比照模式时恢复角向分布曲线选项
        window.ElectronCloud.UI.enableAngularPlotOption();

        // 【新增】退出比照模式时恢复实验面板的原子选择
        const atomSubmenu = document.getElementById('atom-submenu');
        const atomSelectMain = document.getElementById('atom-select');
        if (atomSubmenu) {
            atomSubmenu.classList.remove('disabled');
            atomSubmenu.style.pointerEvents = '';
            atomSubmenu.style.opacity = '';
            atomSubmenu.title = '';
        }
        if (atomSelectMain) {
            atomSelectMain.disabled = false;
        }

        if (label) {
            label.textContent = '选择轨道';
            label.style.display = '';
        }
        if (controlPanel) {
            controlPanel.classList.remove('multiselect-active');
            controlPanel.classList.remove('compare-active');
        }
        if (multiselectControls) multiselectControls.classList.remove('visible');

        // 【新增】隐藏比照模式专用选择器，恢复原有列表
        window.ElectronCloud.UI.toggleCompareOrbitalSelectors(false);
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = false; // 关闭多选模式
            ui.orbitalSelect.size = 3; // 压缩后的高度
        }

        // 清除强制样式类（包括force-unselected）
        if (ui.orbitalSelect) {
            const options = ui.orbitalSelect.querySelectorAll('option');
            options.forEach(option => {
                option.classList.remove('force-selected', 'force-unselected');
                option.style.backgroundColor = '';
                option.style.color = '';
            });
        }

        // 【修复】退出比照模式时，处理选中状态
        if (ui.orbitalSelect) {
            const state = window.ElectronCloud.state;
            const selectedCount = ui.orbitalSelect.selectedOptions.length;
            if (selectedCount > 1) {
                // 多个选中时，保留第一个选中的选项
                const firstSelected = ui.orbitalSelect.selectedOptions[0];
                Array.from(ui.orbitalSelect.options).forEach(option => {
                    option.selected = (option === firstSelected);
                });
            } else if (selectedCount === 0) {
                // 【修改】没有选中时保持空选择状态，与启动时行为一致
                // 不再默认选中第一个选项
            }
            // 同步state.currentOrbitals
            const selected = Array.from(ui.orbitalSelect.selectedOptions).map(o => o.value);
            state.currentOrbitals = selected;

            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);
        }

        // 更新清空按钮状态
        window.ElectronCloud.UI.updateClearAllSelectionsState();
    }
    window.ElectronCloud.UI.updateAngular3DToggleState();
};

// 设置多选模式交互
window.ElectronCloud.UI.setupMultiselectMode = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    if (!ui.orbitalSelect) return;

    // 标记是否已经设置过事件监听器，避免重复设置
    if (ui.orbitalSelect.multiselectHandlerSet) return;
    ui.orbitalSelect.multiselectHandlerSet = true;

    // 存储滚动位置
    let savedScrollTop = 0;

    // 【新增】阻止浏览器默认的focus选择行为
    ui.orbitalSelect.addEventListener('focus', function (event) {
        // 在多选或比照模式下，focus时不应该自动选择任何选项
        if ((ui.multiselectToggle && ui.multiselectToggle.checked) ||
            (ui.compareToggle && ui.compareToggle.checked)) {
            savedScrollTop = ui.orbitalSelect.scrollTop;
        }
    });

    // 多选模式和比照模式的点击处理
    ui.orbitalSelect.addEventListener('mousedown', function (event) {
        const option = event.target;
        if (option.tagName === 'OPTION') {
            // 只有渲染完成后才能切换轨道可见性
            if (state.renderingCompleted && ui.compareToggle && ui.compareToggle.checked) {
                event.preventDefault();
                window.ElectronCloud.Orbital.toggleOrbitalVisibility(option.value);
                return;
            }

            // 只在多选模式或比照模式激活时处理（渲染期间）
            if ((!ui.multiselectToggle || !ui.multiselectToggle.checked) &&
                (!ui.compareToggle || !ui.compareToggle.checked)) return;

            // 保存当前滚动位置
            savedScrollTop = ui.orbitalSelect.scrollTop;

            // 比照模式下限制选择数量
            if (ui.compareToggle && ui.compareToggle.checked && !option.selected && ui.orbitalSelect.selectedOptions.length >= 3) {
                event.preventDefault();
                // 移除最早选中的一个，允许新选择
                const firstSelected = ui.orbitalSelect.selectedOptions[0];
                if (firstSelected) firstSelected.selected = false;
                option.selected = true;
                // alert('比照模式下最多选择三个轨道'); // 不再弹窗，改为自动替换
                return;
            }

            // 阻止默认行为，手动处理选择
            event.preventDefault();

            // 切换选中状态
            option.selected = !option.selected;

            // 立即触发change事件
            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);

            // 立即恢复滚动位置
            setTimeout(() => {
                ui.orbitalSelect.scrollTop = savedScrollTop;
            }, 0);

            // 添加点击动画效果
            option.style.animation = 'selectPulse 0.3s ease-out';
            setTimeout(() => {
                option.style.animation = '';
            }, 300);
        }
    });

    // 【新增】阻止键盘导航自动选择（在多选/比照模式下）
    ui.orbitalSelect.addEventListener('keydown', function (event) {
        if ((ui.multiselectToggle && ui.multiselectToggle.checked) ||
            (ui.compareToggle && ui.compareToggle.checked)) {
            // 阻止方向键等导航键的默认行为（会自动选择选项）
            if (['ArrowUp', 'ArrowDown', 'Home', 'End', 'PageUp', 'PageDown'].includes(event.key)) {
                event.preventDefault();
            }
        }
    });
};

// 更新选择计数
window.ElectronCloud.UI.updateSelectionCount = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    if (!ui.orbitalSelect) return;

    const selectedCount = ui.orbitalSelect.selectedOptions.length;
    const label = document.getElementById('orbital-select-label');
    const countText = document.getElementById('selection-count-text');

    // 更新计数文本
    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        if (countText) countText.textContent = `${selectedCount}`;
        if (label) label.style.display = 'none'; // 隐藏原label
    } else if (ui.compareToggle && ui.compareToggle.checked) {
        if (countText) countText.textContent = `${selectedCount}/3`;
        if (label) label.style.display = 'none'; // 隐藏原label
    } else {
        if (label) label.style.display = ''; // 显示原label
    }

    // 强制刷新选择框样式
    window.ElectronCloud.UI.refreshSelectStyles();
};

// 刷新选择框样式
window.ElectronCloud.UI.refreshSelectStyles = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;

    if (!ui.orbitalSelect) return;

    // 强制重新应用样式
    const options = ui.orbitalSelect.querySelectorAll('option');
    options.forEach(option => {
        // 清除所有颜色类，但保留隐藏/可见类（如果渲染完成）
        option.classList.remove('force-selected', 'force-unselected', 'compare-color-0', 'compare-color-1', 'compare-color-2');

        // 如果渲染未完成，清除隐藏/可见类
        if (!state.renderingCompleted) {
            option.classList.remove('orbital-hidden', 'orbital-visible');
            option.title = '';
        }

        if (option.selected) {
            option.classList.add('force-selected');
            if (ui.compareToggle && ui.compareToggle.checked) {
                // 应用颜色方案
                const selectedOptions = Array.from(ui.orbitalSelect.selectedOptions);
                const colorIndex = selectedOptions.indexOf(option);
                if (colorIndex < constants.compareColors.length) {
                    const color = constants.compareColors[colorIndex].value;
                    option.style.backgroundColor = `rgba(${color.map(v => v * 255).join(',')}, 0.8)`;
                    option.style.color = '#000';
                    // 添加颜色类用于CSS样式
                    option.classList.add(`compare-color-${colorIndex}`);
                }
            } else {
                // 清除比照模式颜色
                option.style.backgroundColor = '';
                option.style.color = '';
            }
        } else {
            // 【关键修复】为未选中的选项添加force-unselected类
            // 这样可以覆盖浏览器的:checked伪类幽灵效果
            option.classList.add('force-unselected');
            option.style.backgroundColor = '';
            option.style.color = '';
        }
    });
};

// 清除所有选择
window.ElectronCloud.UI.clearAllSelections = function () {
    const ui = window.ElectronCloud.ui;

    if (!ui.orbitalSelect || !ui.clearAllSelectionsBtn) return;

    // 使用公共函数清除所有选中状态
    window.ElectronCloud.UI.resetOrbitalSelections();

    // 添加清除动画效果
    ui.clearAllSelectionsBtn.style.transform = 'scale(0.95)';
    setTimeout(() => {
        ui.clearAllSelectionsBtn.style.transform = 'scale(1)';
    }, 150);

    // 触发change事件更新选择状态（保留用于更新UI）
    const changeEvent = new Event('change', { bubbles: true });
    ui.orbitalSelect.dispatchEvent(changeEvent);

    // 更新选择计数显示
    window.ElectronCloud.UI.updateSelectionCount();

    // 更新模式切换栏状态（清空后恢复可用）
    window.ElectronCloud.UI.updateModeSwitcherState();

    // 强制刷新样式
    window.ElectronCloud.UI.refreshSelectStyles();
};

// 专门禁用3D角向形状开关（对比模式下使用）
window.ElectronCloud.UI.disableAngular3DToggle = function () {
    const ui = window.ElectronCloud.ui;

    // 确定当前的模式
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    const isMultiselectMode = ui.multiselectToggle && ui.multiselectToggle.checked;
    const modeText = isCompareMode ? '对比模式下不可用' : isMultiselectMode ? '多选模式下不可用' : '当前模式下不可用';

    // 只禁用3D角向形状开关
    if (ui.angular3dToggle) {
        ui.angular3dToggle.disabled = true;
        ui.angular3dToggle.checked = false;
        const angular3dLabel = document.querySelector('label[for="angular-3d-toggle"]');
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '0.5';
            angular3dLabel.style.cursor = 'not-allowed';
            angular3dLabel.title = modeText;
        }
        // 触发事件以清除角向叠加
        const changeEvent = new Event('change', { bubbles: true });
        ui.angular3dToggle.dispatchEvent(changeEvent);
    }
};

// 专门启用3D角向形状开关（退出对比模式时使用）
window.ElectronCloud.UI.enableAngular3DToggle = function () {
    const ui = window.ElectronCloud.ui;

    // 启用3D角向形状开关，但需要检查其他条件
    if (ui.angular3dToggle) {
        const angular3dLabel = document.querySelector('label[for="angular-3d-toggle"]');
        if (angular3dLabel) {
            angular3dLabel.style.opacity = '1';
            angular3dLabel.style.cursor = 'pointer';
            angular3dLabel.title = '';
        }
        // 调用现有的状态更新函数来正确设置可用性
        window.ElectronCloud.UI.updateAngular3DToggleState();
    }
};

// 禁用角向分布曲线选项（杂化模式/多选模式下使用）
window.ElectronCloud.UI.disableAngularPlotOption = function () {
    const ui = window.ElectronCloud.ui;

    if (ui.plotTypeSelect) {
        // 如果当前选中的是角向分布，切换到径向分布
        if (ui.plotTypeSelect.value === 'angular') {
            ui.plotTypeSelect.value = 'radial';
            // 触发change事件更新图表
            const changeEvent = new Event('change', { bubbles: true });
            ui.plotTypeSelect.dispatchEvent(changeEvent);
        }

        // 禁用整个 select 控件
        ui.plotTypeSelect.disabled = true;

        // 禁用自定义下拉框的视觉效果
        const customContainer = ui.plotTypeSelect.nextElementSibling;
        if (customContainer && customContainer.classList.contains('custom-select-container')) {
            customContainer.style.opacity = '0.5';
            customContainer.style.pointerEvents = 'none';
        }
    }
};

// 启用角向分布曲线选项（退出杂化/多选模式时使用）
window.ElectronCloud.UI.enableAngularPlotOption = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    // 【关键】杂化模式下不启用，防止被其他调用覆盖禁用状态
    if (state.isHybridMode === true) {
        return;
    }

    if (ui.plotTypeSelect) {
        // 启用整个 select 控件
        ui.plotTypeSelect.disabled = false;

        // 恢复自定义下拉框的视觉效果
        const customContainer = ui.plotTypeSelect.nextElementSibling;
        if (customContainer && customContainer.classList.contains('custom-select-container')) {
            customContainer.style.opacity = '';
            customContainer.style.pointerEvents = '';
        }
    }
};

// 更新清空所有选择按钮的状态
window.ElectronCloud.UI.updateClearAllSelectionsState = function () {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    if (!ui.clearAllSelectionsBtn) return;

    // 在对比模式下，如果已经启动计算，则禁用清空按钮
    const isCompareMode = ui.compareToggle && ui.compareToggle.checked;
    const hasStartedDrawing = state.pointCount > 0 || state.isDrawing;

    if (isCompareMode && hasStartedDrawing) {
        ui.clearAllSelectionsBtn.disabled = true;
        ui.clearAllSelectionsBtn.style.opacity = '0.5';
        ui.clearAllSelectionsBtn.style.cursor = 'not-allowed';
        ui.clearAllSelectionsBtn.style.backgroundColor = '';
        ui.clearAllSelectionsBtn.style.color = '';
        ui.clearAllSelectionsBtn.title = '对比模式下启动计算后需要点击重置才能清空选择';
    } else {
        ui.clearAllSelectionsBtn.disabled = false;
        ui.clearAllSelectionsBtn.style.opacity = '1';
        ui.clearAllSelectionsBtn.style.cursor = 'pointer';
        // 可用时标红
        ui.clearAllSelectionsBtn.style.backgroundColor = 'rgba(255, 80, 80, 0.8)';
        ui.clearAllSelectionsBtn.style.color = '#fff';
        ui.clearAllSelectionsBtn.title = '';
    }
};

// 更新全屏按钮状态
window.ElectronCloud.UI.updateFullscreenBtnState = function () {
    const fullscreenBtn = document.getElementById('app-fullscreen-btn');

    if (!fullscreenBtn) return;

    // 允许在任意时候切换全屏，不再在渲染中禁用
    fullscreenBtn.disabled = false;
    fullscreenBtn.title = document.fullscreenElement ? '退出全屏' : '全屏模式';
};

// 启用手势按钮
window.ElectronCloud.UI.enableGestureButton = function () {
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.disabled = false;
        btn.style.opacity = '1';
        btn.style.cursor = 'pointer';
        btn.title = '手势控制';
    }
};

// 禁用手势按钮
window.ElectronCloud.UI.disableGestureButton = function () {
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.disabled = true;
        btn.style.opacity = '0.5';
        btn.style.cursor = 'not-allowed';
        btn.title = '手势控制 (渲染完成后可用)';

        // 如果正在运行手势，停止它
        if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.stop) {
            window.ElectronCloud.Gesture.stop();
        }
    }
};

// 手势按钮点击处理 - 支持开始和停止切换
window.ElectronCloud.UI.onGestureButtonClick = function () {
    // 检查手势识别是否正在运行
    if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.isRunning && window.ElectronCloud.Gesture.isRunning()) {
        // 正在运行，则停止
        window.ElectronCloud.Gesture.stop();

        // 退出全屏
        if (document.fullscreenElement) {
            document.exitFullscreen().catch(err => {
                console.error('退出全屏失败:', err);
            });
        }

        // 更新按钮状态
        const btn = document.getElementById('gesture-control-btn');
        if (btn) {
            btn.classList.remove('gesture-active');
            btn.title = '手势控制';
        }

        return;
    }

    // 未运行，则启动
    // 1. 进入全屏
    if (!document.fullscreenElement) {
        document.documentElement.requestFullscreen().catch(err => {
            console.error(`Error attempting to enable full-screen mode: ${err.message} (${err.name})`);
        });
    }

    // 2. 收起所有面板
    const controlPanel = document.getElementById('control-panel');
    const dataPanel = document.getElementById('data-panel');
    const viewPanel = document.getElementById('view-panel');
    const experimentalPanel = document.getElementById('experimental-panel');
    if (controlPanel) controlPanel.classList.add('collapsed');
    if (dataPanel) dataPanel.classList.add('collapsed');
    if (viewPanel) viewPanel.classList.add('collapsed');
    if (experimentalPanel) experimentalPanel.classList.add('collapsed');

    // 3. 清除视图面板的所有设置
    const state = window.ElectronCloud.state;

    // 3.1 解锁视角
    const viewLockBtns = document.querySelectorAll('.view-lock-btn');
    viewLockBtns.forEach(btn => btn.classList.remove('active'));
    window.ElectronCloud.UI.clearLockedAxisRotation(); // 清除自定义旋转处理
    if (state.controls) {
        state.controls.enabled = true;
        // 注意：不设置 enableRotate，因为使用自定义轨迹球旋转
        state.lockedAxis = null;
        // 恢复默认up向量
        if (state.camera) state.camera.up.set(0, 1, 0);
    }

    // 3.2 停止自动旋转（但保持默认速度为10%）
    if (state.autoRotate) {
        state.autoRotate.enabled = false;
        state.autoRotate.speed = 1; // 默认速度为最大值的10%
    }
    const rotationSpeedRange = document.getElementById('rotation-speed-range');
    const rotationSpeedValue = document.getElementById('rotation-speed-value');
    if (rotationSpeedRange) rotationSpeedRange.value = 1;
    if (rotationSpeedValue) rotationSpeedValue.textContent = '1.00';

    // 3.3 重置旋转轴为默认值
    const rotationAxisX = document.getElementById('rotation-axis-x');
    const rotationAxisY = document.getElementById('rotation-axis-y');
    const rotationAxisZ = document.getElementById('rotation-axis-z');
    if (rotationAxisX) rotationAxisX.value = 0;
    if (rotationAxisY) rotationAxisY.value = 1;
    if (rotationAxisZ) rotationAxisZ.value = 0;
    if (state.autoRotate && state.autoRotate.axis) {
        state.autoRotate.axis.set(0, 1, 0);
    }

    // 4. 开启"固定至中心"
    const centerLock = document.getElementById('center-lock');
    if (centerLock && !centerLock.checked) {
        centerLock.checked = true;
        // 触发 change 事件以应用逻辑
        centerLock.dispatchEvent(new Event('change'));
    }

    // 5. 调整视角 (缩放至合适大小)
    if (state.controls) {
        if (state.farthestDistance > 0) {
            const targetDist = state.farthestDistance * 2.5;
            const vec = new THREE.Vector3().copy(state.camera.position).sub(state.controls.target);
            vec.normalize().multiplyScalar(targetDist);
            state.camera.position.copy(state.controls.target).add(vec);
            state.controls.update();
        }
    }

    // 6. 更新按钮状态
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.classList.add('gesture-active');
        btn.title = '点击停止手势控制';
    }

    // 7. 开始手势识别
    if (window.ElectronCloud.Gesture && window.ElectronCloud.Gesture.start) {
        window.ElectronCloud.Gesture.start();
    } else {
        console.error("Gesture module not loaded");
    }
};

// 更新自动旋转按钮状态
// 【修改】允许在渲染启动前点击开关预设状态，渲染后才实际旋转
window.ElectronCloud.UI.updateAutoRotateButtonState = function () {
    const autoRotateToggle = document.getElementById('auto-rotate-toggle');
    const rotationFeatureBox = document.getElementById('rotation-feature-box');

    // 始终允许用户切换（除非视角被锁定，那由互斥逻辑处理）
    // 实际旋转逻辑在scene.js中已有state.points检查，确保渲染后才旋转
    if (autoRotateToggle) {
        autoRotateToggle.disabled = false;
    }
    if (rotationFeatureBox) {
        rotationFeatureBox.style.opacity = '';
        rotationFeatureBox.style.pointerEvents = '';
    }

    // 同时更新录制按钮状态
    window.ElectronCloud.UI.updateRecordButtonState();
};

// 更新录制按钮状态
window.ElectronCloud.UI.updateRecordButtonState = function () {
    const state = window.ElectronCloud.state;
    const recordBtn = document.getElementById('record-rotation-btn');

    if (!recordBtn) return;

    // 录制按钮启用条件：有点云 && 自动旋转已启用
    const canRecord = state.points &&
        state.autoRotate &&
        state.autoRotate.enabled === true;

    // 如果正在录制，保持启用状态
    if (state.isRecordingRotation) {
        recordBtn.disabled = false;
    } else {
        recordBtn.disabled = !canRecord;
    }
};

// 切换录制状态
window.ElectronCloud.UI.toggleRotationRecording = function () {
    const state = window.ElectronCloud.state;

    if (state.isRecordingRotation) {
        // 停止录制
        window.ElectronCloud.UI.stopRotationRecording();
    } else {
        // 开始录制
        window.ElectronCloud.UI.startRotationRecording();
    }
};

// 开始录制旋转视频
window.ElectronCloud.UI.startRotationRecording = function () {
    const state = window.ElectronCloud.state;
    const recordBtn = document.getElementById('record-rotation-btn');

    if (!state.renderer || !state.autoRotate || !state.autoRotate.enabled) {
        alert('请先启用自动旋转');
        return;
    }

    // 检查浏览器支持
    if (!window.MediaRecorder) {
        alert('您的浏览器不支持视频录制功能');
        return;
    }

    // 【全屏画布方案】画布已经是全屏尺寸，直接录制即可
    const canvas = state.renderer.domElement;

    // 尝试获取支持的 MIME 类型
    let mimeType = 'video/webm;codecs=vp9';
    if (!MediaRecorder.isTypeSupported(mimeType)) {
        mimeType = 'video/webm;codecs=vp8';
        if (!MediaRecorder.isTypeSupported(mimeType)) {
            mimeType = 'video/webm';
            if (!MediaRecorder.isTypeSupported(mimeType)) {
                alert('您的浏览器不支持 WebM 视频录制');
                return;
            }
        }
    }

    try {
        // 获取画布流，60fps
        const stream = canvas.captureStream(60);

        // 计算视频比特率 - 基于实际像素分辨率
        // 高分辨率需要更高比特率才能保持质量
        const pixelCount = canvas.width * canvas.height;
        // 基准: 1080p (1920x1080 = 2073600 像素) 用 12Mbps
        // 按像素数线性缩放，最低 8Mbps，最高 30Mbps
        const baseBitrate = 12000000; // 12 Mbps for 1080p
        const basePixels = 1920 * 1080;
        let videoBitrate = Math.round(baseBitrate * (pixelCount / basePixels));
        videoBitrate = Math.max(8000000, Math.min(videoBitrate, 30000000));

        // 创建 MediaRecorder
        state.mediaRecorder = new MediaRecorder(stream, {
            mimeType: mimeType,
            videoBitsPerSecond: videoBitrate
        });

        state.recordedChunks = [];

        state.mediaRecorder.ondataavailable = (e) => {
            if (e.data.size > 0) {
                state.recordedChunks.push(e.data);
            }
        };

        state.mediaRecorder.onstop = () => {
            // 生成视频文件
            const blob = new Blob(state.recordedChunks, { type: mimeType });
            const url = URL.createObjectURL(blob);
            const link = document.createElement('a');
            const timestamp = new Date().toISOString().replace(/[:.]/g, '-');
            link.download = `electron-cloud-rotation-${timestamp}.webm`;
            link.href = url;
            document.body.appendChild(link);
            link.click();
            setTimeout(() => {
                document.body.removeChild(link);
                URL.revokeObjectURL(url);
            }, 100);

            // 提示用户如何使用WebM文件
            setTimeout(() => {
                alert('视频已保存为 WebM 格式\n\n使用方法：\n• 可直接在浏览器中播放\n• 可上传到 B站、YouTube 等平台\n• 如需转为 MP4，可使用在线工具或 FFmpeg');
            }, 200);
        };

        // 【关键】重置累计角度，从0开始计算
        state.autoRotate.totalAngle = 0;
        state.isRecordingRotation = true;

        // 开始录制
        state.mediaRecorder.start(100); // 每100ms收集一次数据

        // 更新按钮状态
        if (recordBtn) {
            recordBtn.textContent = '停止';
            recordBtn.classList.add('recording');
        }

    } catch (err) {
        console.error('录制启动失败:', err);
        alert('录制启动失败: ' + err.message);
    }
};

// 停止录制旋转视频
window.ElectronCloud.UI.stopRotationRecording = function () {
    const state = window.ElectronCloud.state;
    const recordBtn = document.getElementById('record-rotation-btn');

    if (state.mediaRecorder && state.mediaRecorder.state !== 'inactive') {
        state.mediaRecorder.stop();
    }

    state.isRecordingRotation = false;
    state.autoRotate.totalAngle = 0;

    // 【全屏画布方案】不需要恢复尺寸，画布始终是全屏尺寸

    // 更新按钮状态
    if (recordBtn) {
        recordBtn.textContent = '录制';
        recordBtn.classList.remove('recording');
    }
};

// ========== 自定义轨道列表逻辑 ==========

window.ElectronCloud.UI.initCustomOrbitalList = function () {
    const select = document.getElementById('orbital-select');
    const container = document.getElementById('custom-orbital-list');
    if (!select || !container) return;

    // 清空容器
    container.innerHTML = '';

    // 遍历原生 select 的 options 生成自定义项
    Array.from(select.options).forEach((option, index) => {
        const item = document.createElement('div');
        item.className = 'glass-item';
        item.dataset.value = option.value;
        item.dataset.index = index;

        // 文本内容（不再创建multi-indicator）
        const text = document.createElement('span');
        text.textContent = option.text;
        item.appendChild(text);

        // 点击事件
        item.addEventListener('click', (e) => {
            window.ElectronCloud.UI.handleCustomItemClick(option.value, index, e);
        });

        container.appendChild(item);
    });

    // 【关键修复】初始化时显式清除所有选中状态
    // 确保单选模式与多选/比照模式一样，启动时没有默认选中的轨道
    Array.from(select.options).forEach(option => {
        option.selected = false;
    });
    window.ElectronCloud.state.currentOrbitals = [];

    // 初始化视觉状态
    window.ElectronCloud.UI.updateCustomListVisuals();

    // 初始化时根据当前原子过滤
    const currentAtom = window.ElectronCloud.state.currentAtom || 'H';
    if (window.ElectronCloud.UI.updateOrbitalSelector) {
        window.ElectronCloud.UI.updateOrbitalSelector(currentAtom);
    }

    // 监听原生 select 变化
    select.addEventListener('change', () => {
        window.ElectronCloud.UI.updateCustomListVisuals();
    });
};

window.ElectronCloud.UI.handleCustomItemClick = function (value, index, event) {
    const select = document.getElementById('orbital-select');
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;

    // 渲染完成后，比照模式下允许开关轨道显示
    if (state.renderingCompleted && ui.compareToggle && ui.compareToggle.checked) {
        // 调用toggleOrbitalVisibility切换轨道可见性，而不是改变选中状态
        window.ElectronCloud.Orbital.toggleOrbitalVisibility(value);
        // 更新自定义列表的视觉状态
        window.ElectronCloud.UI.updateCustomListVisualsWithVisibility();
        return;
    }

    // 渲染开始后（有点或正在绘制）但未完成时，禁止操作
    const hasStarted = state.pointCount > 0 || state.isDrawing;
    if (hasStarted) {
        // 渲染进行中禁止操作
        return;
    }

    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        select.options[index].selected = !select.options[index].selected;
    } else if (ui.compareToggle && ui.compareToggle.checked) {
        const selectedOptions = Array.from(select.selectedOptions);
        if (select.options[index].selected) {
            select.options[index].selected = false;
        } else {
            // 允许最多选择3个轨道 (对应红绿蓝)
            if (selectedOptions.length >= 3) {
                // 移除最早选中的一个 (FIFO)
                selectedOptions[0].selected = false;
            }
            select.options[index].selected = true;
        }
    } else {
        Array.from(select.options).forEach(opt => opt.selected = false);
        select.options[index].selected = true;
    }

    select.dispatchEvent(new Event('change'));
    window.ElectronCloud.UI.updateCustomListVisuals();
};

window.ElectronCloud.UI.updateCustomListVisuals = function () {
    const select = document.getElementById('orbital-select');
    const container = document.getElementById('custom-orbital-list');
    const ui = window.ElectronCloud.ui;
    if (!select || !container) return;

    const isMulti = ui.multiselectToggle && ui.multiselectToggle.checked;
    const isCompare = ui.compareToggle && ui.compareToggle.checked;
    const items = container.children;
    const selectedOptions = Array.from(select.selectedOptions);

    Array.from(items).forEach((item, index) => {
        const option = select.options[index];
        const isSelected = option.selected;
        // 移除所有状态类，保持统一
        item.classList.remove('active', 'compare-a', 'compare-b', 'compare-c');

        if (isCompare) {
            // 比照模式：使用红绿蓝颜色
            if (isSelected) {
                const selectionIndex = selectedOptions.indexOf(option);
                if (selectionIndex === 0) item.classList.add('compare-a');
                else if (selectionIndex === 1) item.classList.add('compare-b');
                else if (selectionIndex === 2) item.classList.add('compare-c');
                else item.classList.add('active');
            }
        } else {
            // 普通模式和多选模式：使用统一的绿色active样式
            if (isSelected) item.classList.add('active');
        }
    });
};

// 根据轨道可见性更新自定义列表视觉（渲染完成后比照模式专用）
window.ElectronCloud.UI.updateCustomListVisualsWithVisibility = function () {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    const container = document.getElementById('custom-orbital-list');
    const select = document.getElementById('orbital-select');
    if (!container || !select) return;

    // 仅在比照模式且渲染完成后生效
    if (!state.renderingCompleted || !ui.compareToggle || !ui.compareToggle.checked) {
        window.ElectronCloud.UI.updateCustomListVisuals();
        return;
    }

    const items = container.children;
    const constants = window.ElectronCloud.constants;

    Array.from(items).forEach((item, index) => {
        const option = select.options[index];
        const orbitalKey = option.value;
        const isRendered = state.currentOrbitals && state.currentOrbitals.includes(orbitalKey);
        const isVisible = isRendered && state.orbitalVisibility[orbitalKey] !== false;

        // 清除所有状态类
        item.classList.remove('active', 'compare-a', 'compare-b', 'compare-c', 'orbital-hidden');

        if (isRendered) {
            // 使用渲染时分配的颜色
            const colorIndex = state.currentOrbitals.indexOf(orbitalKey);
            if (isVisible) {
                // 可见状态：显示对应颜色
                if (colorIndex === 0) item.classList.add('compare-a');
                else if (colorIndex === 1) item.classList.add('compare-b');
                else if (colorIndex === 2) item.classList.add('compare-c');
                else item.classList.add('active');
            } else {
                // 隐藏状态：添加隐藏样式类
                item.classList.add('orbital-hidden');
                // 同时保留颜色标识，但通过opacity降低显示
                if (colorIndex === 0) item.classList.add('compare-a');
                else if (colorIndex === 1) item.classList.add('compare-b');
                else if (colorIndex === 2) item.classList.add('compare-c');
            }
        }
        // 非渲染的轨道保持无样式（灰色）
    });
};

// ========== 通用自定义下拉框逻辑 ==========

/**
 * 【性能优化】全局事件委托处理所有自定义下拉框的关闭
 * 替代原来每个createCustomSelect都注册一个全局click事件的做法
 */
window.ElectronCloud.UI._handleGlobalClickForCustomSelects = function (e) {
    // 如果点击的不是自定义下拉框内部，关闭所有打开的下拉框
    if (!e.target.closest('.custom-select-container')) {
        document.querySelectorAll('.custom-select-container.open').forEach(el => {
            el.classList.remove('open');
        });
    }
};

window.ElectronCloud.UI.initCustomSelects = function () {
    // 查找所有非 orbital-select 的 select 元素（包括 max-points-select）
    // 【重要】排除比照模式选择器，它们由 initCompareOrbitalSelectors 单独处理
    const selects = document.querySelectorAll('select:not(#orbital-select):not(.compare-orbital-select):not(.compare-atom-select)');
    selects.forEach(select => {
        // 避免重复初始化
        if (select.classList.contains('hidden-native')) return;
        window.ElectronCloud.UI.createCustomSelect(select);
    });
};

window.ElectronCloud.UI.createCustomSelect = function (select) {
    // 隐藏原生 select
    select.classList.add('hidden-native');

    // 创建容器
    const container = document.createElement('div');
    container.className = 'custom-select-container';

    // 创建触发器
    const trigger = document.createElement('div');
    trigger.className = 'custom-select-trigger';
    trigger.textContent = select.options[select.selectedIndex].text;

    // 创建选项列表
    const optionsList = document.createElement('div');
    optionsList.className = 'custom-select-options';

    // 【性能优化】使用 DocumentFragment 批量创建DOM，减少重排重绘
    const fragment = document.createDocumentFragment();
    const options = select.options;
    const optionsLen = options.length;

    for (let index = 0; index < optionsLen; index++) {
        const option = options[index];
        const customOption = document.createElement('div');
        customOption.className = 'custom-option';
        if (option.selected) customOption.classList.add('selected');
        customOption.textContent = option.text;
        customOption.dataset.value = option.value;
        customOption.dataset.index = index; // 存储索引用于事件委托
        fragment.appendChild(customOption);
    }
    optionsList.appendChild(fragment);

    // 【性能优化】使用事件委托处理选项点击，而非每个选项单独监听
    optionsList.addEventListener('click', (e) => {
        const customOption = e.target.closest('.custom-option');
        if (!customOption) return;
        e.stopPropagation();

        const index = parseInt(customOption.dataset.index, 10);
        // 更新原生 select
        select.selectedIndex = index;
        select.dispatchEvent(new Event('change'));

        // 更新 UI
        trigger.textContent = select.options[index].text;
        container.classList.remove('open');

        // 更新选中状态样式
        Array.from(optionsList.children).forEach(child => child.classList.remove('selected'));
        customOption.classList.add('selected');
    });

    // 触发器点击事件
    trigger.addEventListener('click', (e) => {
        e.stopPropagation();
        // 【关键】如果原生 select 被禁用，不展开自定义下拉框
        if (select.disabled) {
            return;
        }
        // 关闭其他打开的下拉框
        document.querySelectorAll('.custom-select-container.open').forEach(el => {
            if (el !== container) el.classList.remove('open');
        });
        container.classList.toggle('open');
    });

    // 组装 DOM
    container.appendChild(trigger);
    container.appendChild(optionsList);

    // 插入到原生 select 后面
    select.parentNode.insertBefore(container, select.nextSibling);

    // 【性能优化】移除此处的全局click事件注册，改由init()中统一事件委托处理
    // 避免每个select都注册一个全局监听器
};

// ==================== 顶部模式切换栏 ====================
window.ElectronCloud.UI.initModeSwitcher = function () {
    const switcher = document.getElementById('mode-switcher');
    if (!switcher) return;

    const items = switcher.querySelectorAll('.mode-switcher-item');

    // 初始化自动隐藏功能
    this.initModeSwitcherAutoHide();

    items.forEach(item => {
        item.addEventListener('click', () => {
            // 只要电子云还在屏幕上（有点或正在绘制），就禁止切换模式
            const state = window.ElectronCloud.state;
            const hasStarted = state.pointCount > 0 || state.isDrawing;
            if (hasStarted) {
                return; // 阻止模式切换
            }

            const mode = item.dataset.mode;

            // 更新UI状态
            items.forEach(i => i.classList.remove('active'));
            item.classList.add('active');

            // 根据模式切换逻辑
            const multiselectToggle = document.getElementById('multiselect-toggle');
            const compareToggle = document.getElementById('compare-toggle');
            const controlPanel = document.getElementById('control-panel');

            switch (mode) {
                case 'single':
                    // 如果从比照模式切换过来，必须触发事件以清理比照模式的UI状态
                    if (compareToggle && compareToggle.checked) {
                        compareToggle.checked = false;
                        compareToggle.dispatchEvent(new Event('change'));
                    }
                    if (multiselectToggle) multiselectToggle.checked = false;
                    // 移除杂化模式的类和状态
                    if (controlPanel) controlPanel.classList.remove('hybrid-active');
                    state.isHybridMode = false;
                    // 切换为数据面板模式
                    window.ElectronCloud.UI.setHybridPanelMode(false);
                    // 触发change事件以执行原有逻辑
                    multiselectToggle?.dispatchEvent(new Event('change'));
                    break;
                case 'multi':
                    // 如果从比照模式切换过来，必须触发事件以清理比照模式的UI状态
                    if (compareToggle && compareToggle.checked) {
                        compareToggle.checked = false;
                        compareToggle.dispatchEvent(new Event('change'));
                    }
                    if (multiselectToggle) multiselectToggle.checked = true;
                    // 移除杂化模式的类和状态
                    if (controlPanel) controlPanel.classList.remove('hybrid-active');
                    state.isHybridMode = false;
                    // 切换为数据面板模式
                    window.ElectronCloud.UI.setHybridPanelMode(false);
                    multiselectToggle?.dispatchEvent(new Event('change'));
                    break;
                case 'compare':
                    if (multiselectToggle) multiselectToggle.checked = false;
                    if (compareToggle) compareToggle.checked = true;
                    // 移除杂化模式的类和状态
                    if (controlPanel) controlPanel.classList.remove('hybrid-active');
                    state.isHybridMode = false;
                    // 切换为数据面板模式
                    window.ElectronCloud.UI.setHybridPanelMode(false);
                    compareToggle?.dispatchEvent(new Event('change'));
                    break;
                case 'hybrid':
                    // 如果从比照模式切换过来，必须触发事件以清理比照模式的UI状态
                    if (compareToggle && compareToggle.checked) {
                        compareToggle.checked = false;
                        compareToggle.dispatchEvent(new Event('change'));
                    }
                    // 杂化模式 - 复用多选模式的样式和禁用规则
                    if (multiselectToggle) multiselectToggle.checked = true;
                    // 添加杂化模式的类（用于区分）
                    if (controlPanel) controlPanel.classList.add('hybrid-active');
                    // 【核心】设置杂化模式状态
                    state.isHybridMode = true;
                    // 切换为杂化面板模式
                    window.ElectronCloud.UI.setHybridPanelMode(true);
                    // 重置杂化采样缓存
                    if (window.ElectronCloud.Sampling && window.ElectronCloud.Sampling.resetHybridCache) {
                        window.ElectronCloud.Sampling.resetHybridCache();
                    }
                    // 【修复】从比照模式切换过来时，恢复相位开关的可用状态
                    // 必须在触发 onMultiselectToggle 之前处理，因为 compareToggle 已经被设为 false
                    {
                        const phaseToggle = document.getElementById('phase-toggle');
                        const phaseBox = document.getElementById('phase-box');
                        if (phaseToggle) {
                            phaseToggle.disabled = false;
                        }
                        if (phaseBox) {
                            phaseBox.classList.remove('disabled');
                        }
                    }
                    // 触发多选模式的逻辑（复用选择界面）
                    multiselectToggle?.dispatchEvent(new Event('change'));
                    // 禁用3D角向形状开关
                    window.ElectronCloud.UI.updateAngular3DToggleState();
                    // 禁用数据面板中的角向分布曲线选项
                    window.ElectronCloud.UI.disableAngularPlotOption();
                    // 确保轨道轮廓可用（因为 updateAngular3DToggleState 在杂化模式下 return 不会调用它）
                    window.ElectronCloud.UI.updateContour3DToggleState();
                    // 更新杂化轨道按钮
                    window.ElectronCloud.UI.updateHybridOrbitalButtons();
                    break;
            }
        });
    });
};

// 模式切换栏自动隐藏功能
window.ElectronCloud.UI.initModeSwitcherAutoHide = function () {
    const switcher = document.getElementById('mode-switcher');
    if (!switcher) return;

    let hideTimeout;
    let isFullscreen = false;
    let mouseInTopArea = false;

    // 检查是否为全屏模式
    const checkFullscreen = () => {
        return document.fullscreenElement ||
            document.webkitFullscreenElement ||
            document.mozFullScreenElement ||
            document.msFullscreenElement;
    };

    // 隐藏模式切换栏
    const hideSwitcher = () => {
        switcher.classList.add('hidden');
    };

    // 显示模式切换栏
    const showSwitcher = () => {
        switcher.classList.remove('hidden');
        if (isFullscreen) {
            // 在全屏模式下，重新启动自动隐藏计时器
            resetAutoHideTimer();
        }
    };

    // 重置自动隐藏计时器
    const resetAutoHideTimer = () => {
        if (hideTimeout) {
            clearTimeout(hideTimeout);
        }
        if (isFullscreen) {
            // 使用固定的3秒延迟
            const delay = 3000;
            hideTimeout = setTimeout(hideSwitcher, delay);
        }
    };

    // 监听鼠标移动事件（只在屏幕上缘）
    const handleMouseMove = (e) => {
        if (!isFullscreen) return;

        const wasInTopArea = mouseInTopArea;
        mouseInTopArea = e.clientY <= 20;

        // 检查鼠标是否在屏幕上缘（前20px区域）
        if (mouseInTopArea) {
            // 只有当切换栏是隐藏状态时才显示
            if (switcher.classList.contains('hidden')) {
                showSwitcher();
            }
        } else if (wasInTopArea && !mouseInTopArea) {
            // 鼠标离开顶部区域，重新开始自动隐藏计时器
            resetAutoHideTimer();
        }
    };

    // 监听全屏状态变化
    const handleFullscreenChange = () => {
        const wasFullscreen = isFullscreen;
        isFullscreen = checkFullscreen();

        if (isFullscreen && !wasFullscreen) {
            // 进入全屏模式
            mouseInTopArea = false; // 重置鼠标位置状态
            resetAutoHideTimer();
        } else if (!isFullscreen && wasFullscreen) {
            // 退出全屏模式
            mouseInTopArea = false; // 重置鼠标位置状态
            showSwitcher();
            if (hideTimeout) {
                clearTimeout(hideTimeout);
            }
        }
    };

    // 绑定事件监听器
    document.addEventListener('fullscreenchange', handleFullscreenChange);
    document.addEventListener('webkitfullscreenchange', handleFullscreenChange);
    document.addEventListener('mozfullscreenchange', handleFullscreenChange);
    document.addEventListener('MSFullscreenChange', handleFullscreenChange);
    document.addEventListener('mousemove', handleMouseMove);

    // 初始检查全屏状态
    handleFullscreenChange();
};

// 更新模式切换栏的可用状态（继承轨道点选逻辑）
window.ElectronCloud.UI.updateModeSwitcherState = function () {
    const switcher = document.getElementById('mode-switcher');
    if (!switcher) return;

    const state = window.ElectronCloud.state;
    // 只要电子云还在屏幕上（有点或正在绘制），就禁用
    const hasStarted = state.pointCount > 0 || state.isDrawing;

    const items = switcher.querySelectorAll('.mode-switcher-item');
    items.forEach(item => {
        item.style.pointerEvents = hasStarted ? 'none' : 'auto';
        item.style.opacity = hasStarted ? '0.5' : '1';
        item.title = hasStarted ? '电子云绘制中，暂时无法切换模式' : '';
    });
};

// ==================== 杂化轨道控制面板 ====================

/**
 * 初始化杂化轨道控制面板
 */
window.ElectronCloud.UI.initHybridPanel = function () {
    const state = window.ElectronCloud.state;

    // 初始化杂化轨道相关状态
    state.hybridOrbitalPointsMap = {}; // 记录每个杂化轨道的点索引
    state.hybridOrbitalVisibility = {}; // 记录每个杂化轨道的可见性

    // 监听轨道选择变化，动态更新按钮
    const select = document.getElementById('orbital-select');
    if (select) {
        select.addEventListener('change', () => {
            if (state.isHybridMode) {
                window.ElectronCloud.UI.updateHybridOrbitalButtons();
            }
        });
    }
};

/**
 * 设置杂化面板模式（切换数据面板/杂化面板显示）
 * @param {boolean} isHybridMode - 是否为杂化模式
 */
window.ElectronCloud.UI.setHybridPanelMode = function (isHybridMode) {
    const dataPanel = document.getElementById('data-panel');
    const dataModeContent = document.getElementById('data-mode-content');
    const hybridModeContent = document.getElementById('hybrid-mode-content');
    const panelTitle = document.getElementById('data-panel-title');
    const panelTab = document.getElementById('data-panel-tab');

    if (!dataPanel || !dataModeContent || !hybridModeContent) return;

    if (isHybridMode) {
        // 切换到杂化模式
        dataModeContent.style.display = 'none';
        hybridModeContent.style.display = 'block';
        dataPanel.classList.add('hybrid-mode');
        if (panelTitle) panelTitle.textContent = '杂化';
        if (panelTab) panelTab.textContent = '杂化 ▸';
    } else {
        // 切换到数据模式
        dataModeContent.style.display = 'block';
        hybridModeContent.style.display = 'none';
        dataPanel.classList.remove('hybrid-mode');
        if (panelTitle) panelTitle.textContent = '数据';
        if (panelTab) panelTab.textContent = '数据 ▸';
    }
};

/**
 * 更新杂化轨道按钮（根据选中轨道数量动态生成）
 */
window.ElectronCloud.UI.updateHybridOrbitalButtons = function () {
    const state = window.ElectronCloud.state;
    const container = document.getElementById('hybrid-orbital-controls');
    const select = document.getElementById('orbital-select');

    if (!container || !select || !state.isHybridMode) return;

    // 获取选中的轨道数量
    const selectedOrbitals = Array.from(select.selectedOptions).map(opt => opt.value);
    const numOrbitals = selectedOrbitals.length;

    // 清空容器
    container.innerHTML = '';

    // 重置可见性状态
    state.hybridOrbitalVisibility = {};

    // 为每个杂化轨道创建按钮
    for (let i = 0; i < numOrbitals; i++) {
        const btn = document.createElement('button');
        btn.className = 'hybrid-orbital-btn disabled';
        btn.textContent = `杂化轨道${i + 1}`;
        btn.dataset.index = i;

        // 初始化可见性为true
        state.hybridOrbitalVisibility[i] = true;

        // 渲染未完成时按钮不可点击
        btn.addEventListener('click', () => {
            if (!state.renderingCompleted) return;
            window.ElectronCloud.UI.toggleHybridOrbitalVisibility(i);
        });

        container.appendChild(btn);
    }
};

/**
 * 渲染完成后启用杂化轨道按钮
 */
window.ElectronCloud.UI.enableHybridOrbitalButtons = function () {
    const container = document.getElementById('hybrid-orbital-controls');
    if (!container) return;

    const buttons = container.querySelectorAll('.hybrid-orbital-btn');
    buttons.forEach(btn => {
        btn.classList.remove('disabled');
        btn.classList.add('clickable');
    });
};

/**
 * 切换杂化轨道的可见性
 * @param {number} orbitalIndex - 杂化轨道索引
 */
window.ElectronCloud.UI.toggleHybridOrbitalVisibility = function (orbitalIndex) {
    const state = window.ElectronCloud.state;
    const container = document.getElementById('hybrid-orbital-controls');

    if (!state.points || !state.pointOrbitalIndices) return;

    // 切换可见性状态
    state.hybridOrbitalVisibility[orbitalIndex] = !state.hybridOrbitalVisibility[orbitalIndex];
    const isVisible = state.hybridOrbitalVisibility[orbitalIndex];

    // 更新按钮样式
    if (container) {
        const btn = container.querySelector(`.hybrid-orbital-btn[data-index="${orbitalIndex}"]`);
        if (btn) {
            if (isVisible) {
                btn.classList.remove('hidden-orbital');
            } else {
                btn.classList.add('hidden-orbital');
            }
        }
    }

    // 【关键修复】使用位置隐藏方式（与比照模式保持一致）
    // 位置隐藏不会被其他颜色效果覆盖，更加稳定可靠
    const positions = state.points.geometry.attributes.position.array;

    // 如果还没有备份原始位置，先备份
    if (!state.originalPositions) {
        state.originalPositions = new Float32Array(positions);
    }

    for (let i = 0; i < state.pointCount; i++) {
        const pointHybridIdx = state.pointOrbitalIndices[i];

        // 只处理属于当前切换轨道的点
        if (pointHybridIdx !== orbitalIndex) continue;

        const i3 = i * 3;
        if (isVisible) {
            // 显示：恢复原始位置
            positions[i3] = state.originalPositions[i3];
            positions[i3 + 1] = state.originalPositions[i3 + 1];
            positions[i3 + 2] = state.originalPositions[i3 + 2];
        } else {
            // 隐藏：移动到视野外
            positions[i3] = 10000;
            positions[i3 + 1] = 10000;
            positions[i3 + 2] = 10000;
        }
    }

    state.points.geometry.attributes.position.needsUpdate = true;

    // 【同步更新】如果轨道轮廓开关开启，需要更新等值面
    const hybridContourToggle = document.getElementById('hybrid-contour-toggle');
    if (hybridContourToggle && hybridContourToggle.checked) {
        // 触发等值面重新渲染
        hybridContourToggle.dispatchEvent(new Event('change'));
    }
};

/**
 * 重置杂化面板状态
 */
window.ElectronCloud.UI.resetHybridPanel = function () {
    const state = window.ElectronCloud.state;
    const container = document.getElementById('hybrid-orbital-controls');

    // 重置状态
    state.hybridOrbitalPointsMap = {};
    state.hybridOrbitalVisibility = {};

    // 清空按钮容器
    if (container) {
        container.innerHTML = '';
    }

    // 重新生成按钮（如果在杂化模式）
    if (state.isHybridMode) {
        window.ElectronCloud.UI.updateHybridOrbitalButtons();
    }
};

/**
 * 杂化模式轨道轮廓开关处理
 */
window.ElectronCloud.UI.onHybridContourToggle = function (e) {
    const state = window.ElectronCloud.state;
    const checked = e.target.checked;

    if (!state.isHybridMode) return;

    // 清除现有轮廓
    if (state.hybridContourOverlays) {
        for (const overlay of state.hybridContourOverlays) {
            if (!overlay) continue; // 【关键修复】跳过空对象
            state.scene.remove(overlay);
            overlay.traverse((child) => {
                if (child.geometry) child.geometry.dispose();
                if (child.material) child.material.dispose();
            });
        }
        state.hybridContourOverlays = [];
    }

    if (checked) {
        // 创建独立的杂化轨道轮廓
        if (window.ElectronCloud.Visualization && window.ElectronCloud.Visualization.createHybridContourOverlays) {
            state.hybridContourOverlays = window.ElectronCloud.Visualization.createHybridContourOverlays();

            // 添加到场景并同步旋转
            for (const overlay of state.hybridContourOverlays) {
                // 【关键修复】跳过空对象（隐藏的轨道）
                if (!overlay) continue;

                if (state.points) {
                    overlay.rotation.copy(state.points.rotation);
                    overlay.updateMatrix();
                }
                state.scene.add(overlay);
            }
        }
    }
};

/**
 * 根据原子类型更新轨道选择列表的可用性
 * @param {string} atomSymbol - 原子符号 (H, He, Li, ...)
 */
window.ElectronCloud.UI.updateOrbitalSelector = function (atomSymbol) {
    const container = document.getElementById('custom-orbital-list');
    if (!container) return;

    const atomData = window.SlaterBasis ? window.SlaterBasis[atomSymbol] : null;
    const items = container.children;
    const select = document.getElementById('orbital-select');

    // 遍历所有选项
    Array.from(items).forEach((item, index) => {
        const val = item.dataset.value; // e.g. "2px"
        // 提取主壳层键 (e.g. "2px" -> "2p")
        const match = val.match(/^(\d+[spdfg])/);
        const key = match ? match[1] : val;

        let enabled = true;

        if (atomSymbol === 'H') {
            // 氢原子：所有轨道可用（理论上无穷多，这里指列表中的所有）
            enabled = true;
        } else {
            // 其他原子：仅显示 SlaterBasis 中定义的轨道
            // 通常仅包含占据轨道和低激发态
            if (atomData && atomData.orbitals && atomData.orbitals[key]) {
                enabled = true;
            } else {
                enabled = false;
            }
        }

        if (enabled) {
            item.classList.remove('disabled');
            item.style.opacity = '1';
            item.style.pointerEvents = 'auto';
            item.title = '';
        } else {
            item.classList.add('disabled');
            item.style.opacity = '0.3';
            item.style.pointerEvents = 'none';
            item.title = `该原子(${atomSymbol})暂无此轨道(${key})数据`;

            // 如果该不可用项被选中，则取消选中
            const option = select.options[index];
            if (option.selected) {
                option.selected = false;
                // 触发取消选中的视觉更新
                item.classList.remove('active');
                item.classList.remove('compare-a');
                item.classList.remove('compare-b');
                item.classList.remove('compare-c');
            }
        }
    });
};

// ========== 比照模式专用选择器逻辑 ==========

/**
 * 初始化比照模式的三行轨道/原子选择器
 */
window.ElectronCloud.UI.initCompareOrbitalSelectors = function () {
    const state = window.ElectronCloud.state;
    const mainOrbitalSelect = document.getElementById('orbital-select');

    // 获取所有轨道选项
    const orbitalOptions = Array.from(mainOrbitalSelect.options);

    // 填充每个比照行的轨道下拉框
    const orbitalSelects = document.querySelectorAll('.compare-orbital-select');
    orbitalSelects.forEach(select => {
        // 跳过已初始化的
        if (select.classList.contains('hidden-native')) return;

        // 清空现有选项（保留"无"选项）
        select.innerHTML = '<option value="" selected>无</option>';

        // 添加所有轨道选项
        orbitalOptions.forEach(opt => {
            const newOpt = document.createElement('option');
            newOpt.value = opt.value;
            newOpt.textContent = opt.text;
            select.appendChild(newOpt);
        });

        // 添加change事件（在创建自定义下拉框之前）
        select.addEventListener('change', window.ElectronCloud.UI.onCompareOrbitalChange);

        // 【关键】填充完选项后，为这个select创建自定义下拉框
        window.ElectronCloud.UI.createCustomSelect(select);
    });

    // 为原子选择框创建自定义下拉框
    const atomSelects = document.querySelectorAll('.compare-atom-select');
    atomSelects.forEach(select => {
        // 跳过已初始化的
        if (select.classList.contains('hidden-native')) return;

        select.addEventListener('change', window.ElectronCloud.UI.onCompareAtomChange);

        // 【关键】为原子选择器创建自定义下拉框
        window.ElectronCloud.UI.createCustomSelect(select);
    });

    // 为颜色指示器添加点击事件
    const indicators = document.querySelectorAll('.compare-color-indicator');
    indicators.forEach(indicator => {
        indicator.addEventListener('click', window.ElectronCloud.UI.onCompareIndicatorClick);
    });

    // 初始化状态
    window.ElectronCloud.UI.updateCompareSelectorsState();
};

/**
 * 比照模式轨道选择变更处理
 */
window.ElectronCloud.UI.onCompareOrbitalChange = function (event) {
    const state = window.ElectronCloud.state;
    const select = event.target;
    const index = parseInt(select.dataset.index, 10);
    const orbitalValue = select.value;

    // 更新state中的配置
    state.compareMode.slots[index].orbital = orbitalValue;

    // 【新增】轨道→原子约束：更新该行原子可用性
    window.ElectronCloud.UI.updateCompareAtomAvailability(index);

    // 同步到主选择器
    window.ElectronCloud.UI.syncCompareToMainSelect();

    // 更新视觉状态
    window.ElectronCloud.UI.updateCompareSelectorsState();
};

/**
 * 比照模式原子选择变更处理
 */
window.ElectronCloud.UI.onCompareAtomChange = function (event) {
    const state = window.ElectronCloud.state;
    const select = event.target;
    const index = parseInt(select.dataset.index, 10);
    const atomValue = select.value;

    // 更新state中的配置
    state.compareMode.slots[index].atom = atomValue;

    // 【关键修复】同步到主选择器，确保slotConfigs被更新
    window.ElectronCloud.UI.syncCompareToMainSelect();

    // 更新该行的轨道下拉框可用性（根据原子类型过滤轨道）
    window.ElectronCloud.UI.updateCompareOrbitalAvailability(index);
};

/**
 * 根据原子类型更新某一行的轨道可用性
 */
window.ElectronCloud.UI.updateCompareOrbitalAvailability = function (index) {
    const orbitalSelect = document.querySelector(`.compare-orbital-select[data-index="${index}"]`);
    const atomSelect = document.querySelector(`.compare-atom-select[data-index="${index}"]`);
    if (!orbitalSelect || !atomSelect) return;

    const atomSymbol = atomSelect.value;
    const atomData = window.SlaterBasis ? window.SlaterBasis[atomSymbol] : null;

    Array.from(orbitalSelect.options).forEach((option, optIndex) => {
        if (optIndex === 0) return; // 跳过"无"选项

        const val = option.value;
        const match = val.match(/^(\d+[spdfg])/);
        const key = match ? match[1] : val;

        let enabled = true;
        if (atomSymbol === 'H') {
            enabled = true;
        } else if (atomData && atomData.orbitals && atomData.orbitals[key]) {
            enabled = true;
        } else {
            enabled = false;
        }

        option.disabled = !enabled;
        option.style.color = enabled ? '' : 'rgba(255,255,255,0.3)';
    });

    // 【新增】同步更新自定义下拉框UI
    window.ElectronCloud.UI.updateCompareOrbitalCustomSelect(orbitalSelect);

    // 如果当前选中的轨道被禁用，自动切换到"无"
    if (orbitalSelect.selectedIndex > 0 && orbitalSelect.options[orbitalSelect.selectedIndex].disabled) {
        orbitalSelect.value = '';
        orbitalSelect.dispatchEvent(new Event('change'));
    }
};

/**
 * 更新比照模式轨道选择器的自定义下拉框UI
 */
window.ElectronCloud.UI.updateCompareOrbitalCustomSelect = function (orbitalSelect) {
    const container = orbitalSelect.nextElementSibling;
    if (!container || !container.classList.contains('custom-select-container')) return;

    const customOptions = container.querySelector('.custom-select-options');
    if (!customOptions) return;

    Array.from(orbitalSelect.options).forEach((option, index) => {
        const customOption = customOptions.children[index];
        if (!customOption) return;

        if (option.disabled) {
            customOption.classList.add('disabled');
            customOption.style.opacity = '0.3';
            customOption.style.pointerEvents = 'none';
            customOption.title = '该原子不支持此轨道';
        } else {
            customOption.classList.remove('disabled');
            customOption.style.opacity = '1';
            customOption.style.pointerEvents = 'auto';
            customOption.title = '';
        }
    });

    // 更新触发器文本
    if (orbitalSelect.selectedIndex >= 0) {
        const trigger = container.querySelector('.custom-select-trigger');
        if (trigger) {
            trigger.textContent = orbitalSelect.options[orbitalSelect.selectedIndex].text;
        }
    }
};

/**
 * 颜色指示器点击处理（开关轨道显示）
 * 【重构】使用slot索引控制可见性，避免相同轨道同灭同亮
 */
window.ElectronCloud.UI.onCompareIndicatorClick = function (event) {
    const state = window.ElectronCloud.state;
    const indicator = event.currentTarget;
    const uiSlotIndex = parseInt(indicator.dataset.index, 10);

    // 只有在渲染完成后才能开关轨道显示
    if (!state.renderingCompleted) {
        return;
    }

    const slot = state.compareMode.slots[uiSlotIndex];
    if (!slot.orbital) return; // 没有选择轨道，不处理

    // 【关键修复】将UI slot索引映射到activeSlots数组索引（Worker中的orbitalIndex）
    // 因为Worker中的orbitalIndex是在过滤空slot后的activeSlots数组中的位置
    const activeSlots = state.compareMode.activeSlots || [];
    const workerOrbitalIndex = activeSlots.findIndex(s => s.slotIndex === uiSlotIndex);

    if (workerOrbitalIndex < 0) {
        console.log(`无法找到UI slot ${uiSlotIndex} 对应的activeSlot`);
        return;
    }

    // 【重构】使用slot索引控制可见性，而非轨道键
    // 初始化slotVisibility对象（如果不存在）
    if (!state.compareMode.slotVisibility) {
        state.compareMode.slotVisibility = {};
    }

    // 切换该slot的可见性（使用UI slot索引记录状态）
    const currentVisible = state.compareMode.slotVisibility[uiSlotIndex] !== false;  // 默认为可见
    state.compareMode.slotVisibility[uiSlotIndex] = !currentVisible;

    // 调用场景更新函数来同步点的显示（使用Worker的orbitalIndex来匹配点）
    window.ElectronCloud.Orbital.updateCompareSlotVisibility(workerOrbitalIndex, !currentVisible);

    // 更新视觉状态
    window.ElectronCloud.UI.updateCompareSelectorsState();
};

/**
 * 同步比照模式选择到主选择器
 * 【重构】使用slot索引作为标识，支持相同轨道不同原子
 */
window.ElectronCloud.UI.syncCompareToMainSelect = function () {
    const state = window.ElectronCloud.state;
    const mainSelect = document.getElementById('orbital-select');
    if (!mainSelect) return;

    // 清除所有选中
    Array.from(mainSelect.options).forEach(opt => opt.selected = false);

    // 【重构】构建slot配置数组，每个slot独立（即使轨道相同）
    const activeSlots = [];  // [{slotIndex, orbital, atom}]
    const selectedOrbitalKeys = new Set();  // 用于更新主选择器的显示

    state.compareMode.slots.forEach((slot, index) => {
        if (slot.orbital) {
            activeSlots.push({
                slotIndex: index,
                orbital: slot.orbital,
                atom: slot.atom || 'H'
            });
            selectedOrbitalKeys.add(slot.orbital);
        }
    });

    // 在主选择器中选中这些轨道（仅用于UI显示）
    selectedOrbitalKeys.forEach(orbital => {
        const option = mainSelect.querySelector(`option[value="${orbital}"]`);
        if (option) option.selected = true;
    });

    // 【重构】保存activeSlots到state（供采样使用）
    state.compareMode.activeSlots = activeSlots;

    // 构建用于Worker的配置：每个slot独立的轨道+原子
    // 格式：[{orbital: '2s', atom: 'H'}, {orbital: '2s', atom: 'C'}, ...]
    state.compareMode.slotConfigs = activeSlots.map(s => ({
        orbital: s.orbital,
        atom: s.atom,
        slotIndex: s.slotIndex
    }));

    // 更新 currentOrbitals（用于兼容其他代码）
    state.currentOrbitals = activeSlots.map(s => s.orbital);

    // 更新选择计数
    window.ElectronCloud.UI.updateSelectionCount();

    // 触发change事件更新相关状态
    mainSelect.dispatchEvent(new Event('change'));
};

/**
 * 更新比照选择器的视觉状态
 * 【重构】使用slotVisibility而非orbitalVisibility
 */
window.ElectronCloud.UI.updateCompareSelectorsState = function () {
    const state = window.ElectronCloud.state;

    const indicators = document.querySelectorAll('.compare-color-indicator');
    indicators.forEach(indicator => {
        const slotIndex = parseInt(indicator.dataset.index, 10);
        const slot = state.compareMode.slots[slotIndex];

        // 根据是否有轨道选中来设置激活状态
        if (slot.orbital) {
            indicator.classList.add('active');

            // 【重构】根据slotVisibility设置隐藏状态
            const slotVisibility = state.compareMode.slotVisibility || {};
            if (state.renderingCompleted && slotVisibility[slotIndex] === false) {
                indicator.classList.add('orbital-hidden');
            } else {
                indicator.classList.remove('orbital-hidden');
            }
        } else {
            indicator.classList.remove('active');
            indicator.classList.remove('orbital-hidden');
        }
    });

    // 更新下拉框的禁用状态（渲染中禁用）
    const hasStarted = state.pointCount > 0 || state.isDrawing;
    const isRendering = hasStarted && !state.renderingCompleted;

    const orbitalSelects = document.querySelectorAll('.compare-orbital-select');
    const atomSelects = document.querySelectorAll('.compare-atom-select');

    orbitalSelects.forEach(select => {
        select.disabled = isRendering;
    });
    atomSelects.forEach(select => {
        select.disabled = isRendering;
    });
};

/**
 * 重置比照模式选择器
 */
window.ElectronCloud.UI.resetCompareOrbitalSelectors = function () {
    const state = window.ElectronCloud.state;

    // 重置state中的配置
    state.compareMode.slots = [
        { orbital: '', atom: 'H' },
        { orbital: '', atom: 'H' },
        { orbital: '', atom: 'H' }
    ];
    state.compareMode.orbitalAtomMap = {};

    // 重置UI
    const orbitalSelects = document.querySelectorAll('.compare-orbital-select');
    orbitalSelects.forEach(select => {
        select.value = '';
    });

    const atomSelects = document.querySelectorAll('.compare-atom-select');
    atomSelects.forEach(select => {
        select.value = 'H';
    });

    // 更新视觉状态
    window.ElectronCloud.UI.updateCompareSelectorsState();
};

/**
 * 显示/隐藏比照模式选择器
 */
window.ElectronCloud.UI.toggleCompareOrbitalSelectors = function (show) {
    const container = document.getElementById('compare-orbital-selectors');
    const customList = document.getElementById('custom-orbital-list');
    const multiselectLabel = document.getElementById('multiselect-label');

    if (container) {
        container.style.display = show ? 'flex' : 'none';
    }

    if (customList) {
        customList.style.display = show ? 'none' : 'block';
    }

    // 比照模式下隐藏"已选"计数标签
    if (multiselectLabel) {
        multiselectLabel.style.display = show ? 'none' : '';
    }
};

// ==================== 比照模式双向约束 ====================

/**
 * 【比照模式专用】根据选中的轨道更新原子可用性
 * @param {number} index - 比照模式slot索引 (0, 1, 2)
 */
window.ElectronCloud.UI.updateCompareAtomAvailability = function (index) {
    const orbitalSelect = document.querySelector(`.compare-orbital-select[data-index="${index}"]`);
    const atomSelect = document.querySelector(`.compare-atom-select[data-index="${index}"]`);
    if (!orbitalSelect || !atomSelect) return;

    const orbitalValue = orbitalSelect.value;

    // 如果没有选择轨道（"无"），则所有原子都可用
    if (!orbitalValue) {
        Array.from(atomSelect.options).forEach(option => {
            option.disabled = false;
            option.style.color = '';
        });
        // 更新自定义下拉框
        window.ElectronCloud.UI.updateCompareAtomCustomSelect(atomSelect);
        return;
    }

    // 提取轨道类型 (2px -> 2p, 4d_xy -> 4d)
    const match = orbitalValue.match(/^(\d+[spdfg])/);
    const orbitalKey = match ? match[1] : orbitalValue;

    Array.from(atomSelect.options).forEach(option => {
        const atomSymbol = option.value;
        let enabled = true;

        if (atomSymbol === 'H') {
            // 氢原子：所有轨道可用
            enabled = true;
        } else {
            const atomData = window.SlaterBasis ? window.SlaterBasis[atomSymbol] : null;
            if (!atomData || !atomData.orbitals || !atomData.orbitals[orbitalKey]) {
                enabled = false;
            }
        }

        option.disabled = !enabled;
        option.style.color = enabled ? '' : 'rgba(255,255,255,0.3)';
    });

    // 更新自定义下拉框UI
    window.ElectronCloud.UI.updateCompareAtomCustomSelect(atomSelect);

    // 如果当前选中的原子被禁用，自动切换到H
    if (atomSelect.selectedIndex >= 0 && atomSelect.options[atomSelect.selectedIndex].disabled) {
        atomSelect.value = 'H';
        atomSelect.dispatchEvent(new Event('change'));
    }
};

/**
 * 更新比照模式原子选择器的自定义下拉框UI
 */
window.ElectronCloud.UI.updateCompareAtomCustomSelect = function (atomSelect) {
    const container = atomSelect.nextElementSibling;
    if (!container || !container.classList.contains('custom-select-container')) return;

    const customOptions = container.querySelector('.custom-select-options');
    if (!customOptions) return;

    Array.from(atomSelect.options).forEach((option, index) => {
        const customOption = customOptions.children[index];
        if (!customOption) return;

        if (option.disabled) {
            customOption.classList.add('disabled');
            customOption.style.opacity = '0.3';
            customOption.style.pointerEvents = 'none';
            customOption.title = '该原子不支持选中的轨道';
        } else {
            customOption.classList.remove('disabled');
            customOption.style.opacity = '1';
            customOption.style.pointerEvents = 'auto';
            customOption.title = '';
        }
    });

    // 更新触发器文本
    if (atomSelect.selectedIndex >= 0) {
        const trigger = container.querySelector('.custom-select-trigger');
        if (trigger) {
            trigger.textContent = atomSelect.options[atomSelect.selectedIndex].text;
        }
    }
};
