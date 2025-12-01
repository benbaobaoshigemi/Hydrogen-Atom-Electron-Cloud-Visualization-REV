// UI 控件交互模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.UI = {};

// 初始化UI事件监听
window.ElectronCloud.UI.init = function() {
    const ui = window.ElectronCloud.ui;
    
    // 初始化自定义轨道列表
    if (window.ElectronCloud.UI.initCustomOrbitalList) {
        window.ElectronCloud.UI.initCustomOrbitalList();
    }
    
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
            window.ElectronCloud.Orbital.clearDrawing();
            if (window.DataPanel && window.DataPanel.reset) {
                window.DataPanel.reset();
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
    
    // 确保初始状态下自动旋转按钮禁用（渲染开始前）
    window.ElectronCloud.UI.updateAutoRotateButtonState();
    
    // 确保初始状态下坐标系正确隐藏（farthestDistance=0时）
    if (window.ElectronCloud.ui.axesSizeRange) {
        window.ElectronCloud.UI.onAxesSizeChange({ target: window.ElectronCloud.ui.axesSizeRange });
    }
    
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
            console.log('闪烁模式:', e.target.checked, 'bloomPass:', !!state.bloomPass);
            
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
            
            console.log('闪烁模式切换:', previousMode, '->', e.target.value);
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
    
    console.log('权重模式初始化 - toggle:', weightModeToggle, 'box:', weightFeatureBox);
    
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
            
            console.log('权重模式:', e.target.checked ? '开启' : '关闭');
        });
    } else {
        console.error('weight-mode-toggle 元素未找到');
    }
    
    // 权重模式功能框点击切换
    if (weightFeatureBox) {
        weightFeatureBox.addEventListener('click', (e) => {
            const toggle = document.getElementById('weight-mode-toggle');
            console.log('权重功能框被点击, toggle:', toggle);
            if (toggle) {
                toggle.checked = !toggle.checked;
                toggle.dispatchEvent(new Event('change', { bubbles: true }));
            }
        });
    } else {
        console.error('weight-feature-box 元素未找到');
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
            console.log('导出功能框被点击');
            if (window.ElectronCloud.Visualization && typeof window.ElectronCloud.Visualization.exportImage === 'function') {
                window.ElectronCloud.Visualization.exportImage();
            } else {
                console.error('导出函数不可用');
                alert('导出功能暂时不可用，请刷新页面后重试');
            }
        });
        console.log('导出功能框点击事件已绑定');
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
    
    console.log('UI 事件监听器初始化完成');

};

// 锁定摄像机到指定轴视角
window.ElectronCloud.UI.lockCameraToAxis = function(axis) {
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
    
    switch(axis) {
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
    
    console.log(`视角已锁定到${axis}轴`);
};

// 设置锁定轴的自定义旋转处理
window.ElectronCloud.UI.setupLockedAxisRotation = function(axis) {
    const state = window.ElectronCloud.state;
    const canvas = state.renderer.domElement;
    
    // 先清除之前的监听器
    window.ElectronCloud.UI.clearLockedAxisRotation();
    
    console.log('设置锁定轴旋转处理, axis:', axis);
    
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
    const getRotationAxis = function() {
        switch(state.lockedAxis) {
            case 'x': return new THREE.Vector3(1, 0, 0);
            case 'y': return new THREE.Vector3(0, 1, 0);
            case 'z': return new THREE.Vector3(0, 0, 1);
            default: return null;
        }
    };
    
    // 旋转场景对象
    const rotateSceneObjects = function(angle) {
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
    const handleZoom = function(delta) {
        const zoomSpeed = 0.001;
        const factor = 1 - delta * zoomSpeed;
        const minDistance = 5;
        const maxDistance = 500;
        
        const currentDistance = state.camera.position.length();
        const newDistance = Math.max(minDistance, Math.min(maxDistance, currentDistance * factor));
        
        state.camera.position.normalize().multiplyScalar(newDistance);
    };
    
    // 使用 pointer 事件以获得更好的兼容性（包括触摸）
    const onPointerDown = function(e) {
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
    
    const onPointerMove = function(e) {
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
    
    const onPointerUp = function(e) {
        if (isDragging) {
            isDragging = false;
            try {
                canvas.releasePointerCapture(e.pointerId);
            } catch(err) {
                // 忽略释放捕获的错误
            }
        }
    };
    
    // 滚轮缩放（包括触控板双指缩放）
    const onWheel = function(e) {
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
    
    const onTouchStart = function(e) {
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
    
    const onTouchMove = function(e) {
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
    
    const onTouchEnd = function(e) {
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
    
    console.log('已设置锁定轴旋转处理，OrbitControls 已禁用，使用自定义 pointer 事件');
};

// 清除锁定轴的自定义旋转处理
window.ElectronCloud.UI.clearLockedAxisRotation = function() {
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
        console.log('已清除锁定轴旋转处理');
    }
    
    // 重新启用 OrbitControls
    if (state.controls) {
        state.controls.enabled = true;
    }
};

// 自动旋转与锁定视角互斥逻辑
window.ElectronCloud.UI.updateRotationLockMutualExclusion = function() {
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
window.ElectronCloud.UI.calculateAxesScale = function(scaleFactor) {
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
window.ElectronCloud.UI.resetOrbitalSelections = function() {
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
    
    console.log('[UI] 轨道选择已清除');
};

// 重置所有场景对象的旋转状态（公共函数，避免重复代码）
// 仅重置旋转，不清除锁定状态或其他逻辑
window.ElectronCloud.UI.resetSceneObjectsRotation = function() {
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
window.ElectronCloud.UI.onAxesSizeChange = function(event) {
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
window.ElectronCloud.UI.updateAxesSizeRangeState = function() {
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
window.ElectronCloud.UI.onAngularOverlayToggle = function(event) {
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

// 处理轨道选择变化
window.ElectronCloud.UI.onOrbitalSelectChange = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    const opts = Array.from(ui.orbitalSelect.selectedOptions || []).map(o => o.value);
    
    // 在多选或比照模式下，使用选中的选项
    // 在普通模式下，如果没有选中任何选项，默认使用第一个选项
    if ((ui.multiselectToggle && ui.multiselectToggle.checked) || (ui.compareToggle && ui.compareToggle.checked)) {
        state.currentOrbitals = opts.length ? opts : [];
    } else {
        state.currentOrbitals = opts.length ? opts : [ui.orbitalSelect.value || '1s'];
    }
    
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
    
    // 更新多选模式或比照模式的选择计数
    if ((ui.multiselectToggle && ui.multiselectToggle.checked) || (ui.compareToggle && ui.compareToggle.checked)) {
        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();
        
        // 更新清空按钮状态
        window.ElectronCloud.UI.updateClearAllSelectionsState();
    }
};

// 单选模式下更新选中样式
window.ElectronCloud.UI.updateSingleSelectStyle = function() {
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
window.ElectronCloud.UI.updateAngular3DToggleState = function() {
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
};

// 控件辅助函数
window.ElectronCloud.UI.onSpeedChange = function() {
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
window.ElectronCloud.UI.getPointSize = function() {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;
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
    
    // 获取当前选中轨道中最大的预估距离
    let maxDistance = baseDistance; // 默认值
    
    // 优先从 state.currentOrbitals 获取（适用于所有模式）
    if (state.currentOrbitals && state.currentOrbitals.length > 0) {
        for (const orbitalKey of state.currentOrbitals) {
            const dist = distanceTable[orbitalKey];
            if (dist && dist > maxDistance) {
                maxDistance = dist;
            }
        }
    } 
    // 备用：从 UI 选择框获取
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

window.ElectronCloud.UI.onSizeChange = function() {
    const state = window.ElectronCloud.state;
    
    if (state.points && state.points.material) {
        state.points.material.size = window.ElectronCloud.UI.getPointSize();
    }
};

window.ElectronCloud.UI.onOpacityChange = function() {
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

window.ElectronCloud.UI.onCenterLockChange = function() {
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
    
    console.log('[UI] centerLock 状态变更:', locked ? '锁定' : '解锁', ', enablePan:', state.controls.enablePan);
};

// 辅助函数：检查并同步 centerLock 状态
// 用于其他模块恢复状态时调用
window.ElectronCloud.UI.syncCenterLockState = function() {
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
window.ElectronCloud.UI.onPhaseToggleChange = function() {
    const phaseToggle = document.getElementById('phase-toggle');
    const phaseBox = document.getElementById('phase-box');
    
    if (phaseBox && phaseToggle) {
        phaseBox.classList.toggle('active', phaseToggle.checked);
    }
};

// 为所有 mode-toggle-box 添加点击切换逻辑
window.ElectronCloud.UI.setupModeToggleBoxes = function() {
    const boxes = document.querySelectorAll('.mode-toggle-box');
    boxes.forEach((box, index) => {
        box.addEventListener('click', function(e) {
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

window.ElectronCloud.UI.onPauseToggle = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    if (!state.isDrawing) {
        // 继续采样，但不启用实时图表更新
        state.isDrawing = true;
        if (ui.pauseButton) ui.pauseButton.textContent = '暂停';
        window.ElectronCloud.UI.updateAngular3DToggleState();
        window.ElectronCloud.UI.updateAxesSizeRangeState();
        return;
    }
    // 暂停
    state.isDrawing = false;
    if (ui.pauseButton) ui.pauseButton.textContent = '继续';
    window.ElectronCloud.UI.updateAngular3DToggleState();
    window.ElectronCloud.UI.updateAxesSizeRangeState();
};

// 多选模式功能
window.ElectronCloud.UI.onMultiselectToggle = function() {
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
                // 没有选中任何选项时，默认选中第一个选项（1s）
                ui.orbitalSelect.options[0].selected = true;
            }
            // 同步state.currentOrbitals
            const selected = Array.from(ui.orbitalSelect.selectedOptions).map(o => o.value);
            state.currentOrbitals = selected.length ? selected : ['1s'];
            
            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);
        }
    }
    window.ElectronCloud.UI.updateAngular3DToggleState();
};

// 比照模式功能
window.ElectronCloud.UI.onCompareToggle = function() {
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
        
        // 禁用显示开关功能，避免造成bug
        window.ElectronCloud.UI.disableAngular3DToggle();
        
        if (controlPanel) {
            controlPanel.classList.add('multiselect-active');
            controlPanel.classList.add('compare-active');
        }
        if (multiselectControls) multiselectControls.classList.add('visible');
        if (label) label.style.display = 'none';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = true;
            ui.orbitalSelect.size = 5;
        }
        
        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();
        window.ElectronCloud.UI.updateCustomListVisuals();
        
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
        
        // 重新启用显示开关功能
        window.ElectronCloud.UI.enableAngular3DToggle();
        
        if (label) {
            label.textContent = '选择轨道';
            label.style.display = '';
        }
        if (controlPanel) {
            controlPanel.classList.remove('multiselect-active');
            controlPanel.classList.remove('compare-active');
        }
        if (multiselectControls) multiselectControls.classList.remove('visible');
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
                // 没有选中任何选项时，默认选中第一个选项（1s）
                ui.orbitalSelect.options[0].selected = true;
            }
            // 同步state.currentOrbitals
            const selected = Array.from(ui.orbitalSelect.selectedOptions).map(o => o.value);
            state.currentOrbitals = selected.length ? selected : ['1s'];
            
            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);
        }
        
        // 更新清空按钮状态
        window.ElectronCloud.UI.updateClearAllSelectionsState();
    }
    window.ElectronCloud.UI.updateAngular3DToggleState();
};

// 设置多选模式交互
window.ElectronCloud.UI.setupMultiselectMode = function() {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;
    
    if (!ui.orbitalSelect) return;
    
    // 标记是否已经设置过事件监听器，避免重复设置
    if (ui.orbitalSelect.multiselectHandlerSet) return;
    ui.orbitalSelect.multiselectHandlerSet = true;
    
    // 存储滚动位置
    let savedScrollTop = 0;
    
    // 【新增】阻止浏览器默认的focus选择行为
    ui.orbitalSelect.addEventListener('focus', function(event) {
        // 在多选或比照模式下，focus时不应该自动选择任何选项
        if ((ui.multiselectToggle && ui.multiselectToggle.checked) || 
            (ui.compareToggle && ui.compareToggle.checked)) {
            savedScrollTop = ui.orbitalSelect.scrollTop;
        }
    });
    
    // 多选模式和比照模式的点击处理
    ui.orbitalSelect.addEventListener('mousedown', function(event) {
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
    ui.orbitalSelect.addEventListener('keydown', function(event) {
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
window.ElectronCloud.UI.updateSelectionCount = function() {
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
window.ElectronCloud.UI.refreshSelectStyles = function() {
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
window.ElectronCloud.UI.clearAllSelections = function() {
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
    
    // 强制刷新样式
    window.ElectronCloud.UI.refreshSelectStyles();
};

// 专门禁用3D角向形状开关（对比模式下使用）
window.ElectronCloud.UI.disableAngular3DToggle = function() {
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
window.ElectronCloud.UI.enableAngular3DToggle = function() {
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

// 禁用角向分布曲线选项（多选模式下使用）
window.ElectronCloud.UI.disableAngularPlotOption = function() {
    const ui = window.ElectronCloud.ui;
    
    if (ui.plotTypeSelect) {
        // 找到角向分布选项并禁用它
        const angularOption = ui.plotTypeSelect.querySelector('option[value="angular"]');
        if (angularOption) {
            angularOption.disabled = true;
            angularOption.title = '多选模式下不可用';
        }
        
        // 如果当前选中的是角向分布，切换到径向分布
        if (ui.plotTypeSelect.value === 'angular') {
            ui.plotTypeSelect.value = 'radial';
            // 触发change事件更新图表
            const changeEvent = new Event('change', { bubbles: true });
            ui.plotTypeSelect.dispatchEvent(changeEvent);
        }
        
        // 更新标签样式提示
        const plotTypeLabel = document.querySelector('label[for="plot-type-select"]');
        if (plotTypeLabel) {
            plotTypeLabel.title = '多选模式下仅支持径向分布';
        }
    }
};

// 启用角向分布曲线选项（退出多选模式时使用）
window.ElectronCloud.UI.enableAngularPlotOption = function() {
    const ui = window.ElectronCloud.ui;
    
    if (ui.plotTypeSelect) {
        // 启用角向分布选项
        const angularOption = ui.plotTypeSelect.querySelector('option[value="angular"]');
        if (angularOption) {
            angularOption.disabled = false;
            angularOption.title = '';
        }
        
        // 清除标签提示
        const plotTypeLabel = document.querySelector('label[for="plot-type-select"]');
        if (plotTypeLabel) {
            plotTypeLabel.title = '';
        }
    }
};

// 更新清空所有选择按钮的状态
window.ElectronCloud.UI.updateClearAllSelectionsState = function() {
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
window.ElectronCloud.UI.updateFullscreenBtnState = function() {
    const fullscreenBtn = document.getElementById('app-fullscreen-btn');
    
    if (!fullscreenBtn) return;
    
    // 允许在任意时候切换全屏，不再在渲染中禁用
    fullscreenBtn.disabled = false;
    fullscreenBtn.title = document.fullscreenElement ? '退出全屏' : '全屏模式';
};

// 启用手势按钮
window.ElectronCloud.UI.enableGestureButton = function() {
    const btn = document.getElementById('gesture-control-btn');
    if (btn) {
        btn.disabled = false;
        btn.style.opacity = '1';
        btn.style.cursor = 'pointer';
        btn.title = '手势控制';
    }
};

// 禁用手势按钮
window.ElectronCloud.UI.disableGestureButton = function() {
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
window.ElectronCloud.UI.onGestureButtonClick = function() {
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
        
        console.log('[UI] 手势控制已停止');
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
    
    // 3.2 停止自动旋转
    if (state.autoRotate) {
        state.autoRotate.enabled = false;
        state.autoRotate.speed = 0;
    }
    const rotationSpeedRange = document.getElementById('rotation-speed-range');
    const rotationSpeedValue = document.getElementById('rotation-speed-value');
    if (rotationSpeedRange) rotationSpeedRange.value = 0;
    if (rotationSpeedValue) rotationSpeedValue.textContent = '0';
    
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

// 更新自动旋转按钮状态（只有在没有点云时才禁用）
window.ElectronCloud.UI.updateAutoRotateButtonState = function() {
    const state = window.ElectronCloud.state;
    const autoRotateToggle = document.getElementById('auto-rotate-toggle');
    const rotationFeatureBox = document.getElementById('rotation-feature-box');
    
    // 只要有点云就可以启用自动旋转
    const canEnable = state.points !== null && state.points !== undefined;
    
    if (autoRotateToggle) {
        autoRotateToggle.disabled = !canEnable;
    }
    
    if (rotationFeatureBox) {
        if (canEnable) {
            rotationFeatureBox.style.opacity = '';
            rotationFeatureBox.style.pointerEvents = '';
        } else {
            rotationFeatureBox.style.opacity = '0.4';
            rotationFeatureBox.style.pointerEvents = 'none';
        }
    }
    
    // 同时更新录制按钮状态
    window.ElectronCloud.UI.updateRecordButtonState();
};

// 更新录制按钮状态
window.ElectronCloud.UI.updateRecordButtonState = function() {
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
window.ElectronCloud.UI.toggleRotationRecording = function() {
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
window.ElectronCloud.UI.startRotationRecording = function() {
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
    
    console.log('录制画布尺寸 (全屏):', canvas.width, 'x', canvas.height);
    
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
        
        console.log('录制参数: 画布尺寸', canvas.width, 'x', canvas.height, 
                    ', 像素数:', pixelCount, ', 比特率:', (videoBitrate / 1000000).toFixed(1), 'Mbps');
        
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
            
            console.log('录制完成，视频已保存');
            
            // 提示用户如何使用WebM文件
            setTimeout(() => {
                alert('视频已保存为 WebM 格式\n\n使用方法：\n• 可直接在浏览器中播放\n• 可上传到 B站、YouTube 等平台\n• 如需转为 MP4，可使用在线工具或 FFmpeg');
            }, 200);
        };
        
        // 【关键】重置累计角度，从0开始计算
        state.autoRotate.totalAngle = 0;
        state.isRecordingRotation = true;
        
        console.log('开始录制，累计角度已重置为0，目标:', (Math.PI * 2).toFixed(4), '(一周)');
        
        // 开始录制
        state.mediaRecorder.start(100); // 每100ms收集一次数据
        
        // 更新按钮状态
        if (recordBtn) {
            recordBtn.textContent = '停止';
            recordBtn.classList.add('recording');
        }
        
        console.log('开始录制旋转视频，等待旋转一周...');
        
    } catch (err) {
        console.error('录制启动失败:', err);
        alert('录制启动失败: ' + err.message);
    }
};

// 停止录制旋转视频
window.ElectronCloud.UI.stopRotationRecording = function() {
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
    
    console.log('录制已停止');
};

// ========== 自定义轨道列表逻辑 ==========

window.ElectronCloud.UI.initCustomOrbitalList = function() {
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

    // 初始化视觉状态
    window.ElectronCloud.UI.updateCustomListVisuals();
    
    // 监听原生 select 变化
    select.addEventListener('change', () => {
        window.ElectronCloud.UI.updateCustomListVisuals();
    });
};

window.ElectronCloud.UI.handleCustomItemClick = function(value, index, event) {
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

window.ElectronCloud.UI.updateCustomListVisuals = function() {
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
window.ElectronCloud.UI.updateCustomListVisualsWithVisibility = function() {
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

window.ElectronCloud.UI.initCustomSelects = function() {
    // 查找所有非 orbital-select 的 select 元素（包括 max-points-select）
    const selects = document.querySelectorAll('select:not(#orbital-select)');
    selects.forEach(select => {
        // 避免重复初始化
        if (select.classList.contains('hidden-native')) return;
        window.ElectronCloud.UI.createCustomSelect(select);
    });
};

window.ElectronCloud.UI.createCustomSelect = function(select) {
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
    
    // 填充选项
    Array.from(select.options).forEach((option, index) => {
        const customOption = document.createElement('div');
        customOption.className = 'custom-option';
        if (option.selected) customOption.classList.add('selected');
        customOption.textContent = option.text;
        customOption.dataset.value = option.value;
        
        customOption.addEventListener('click', (e) => {
            e.stopPropagation();
            // 更新原生 select
            select.selectedIndex = index;
            select.dispatchEvent(new Event('change'));
            
            // 更新 UI
            trigger.textContent = option.text;
            container.classList.remove('open');
            
            // 更新选中状态样式
            Array.from(optionsList.children).forEach(child => child.classList.remove('selected'));
            customOption.classList.add('selected');
        });
        
        optionsList.appendChild(customOption);
    });
    
    // 触发器点击事件
    trigger.addEventListener('click', (e) => {
        e.stopPropagation();
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
    
    // 点击外部关闭
    document.addEventListener('click', (e) => {
        if (!container.contains(e.target)) {
            container.classList.remove('open');
        }
    });
};