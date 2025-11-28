// UI 控件交互模块
window.ElectronCloud = window.ElectronCloud || {};
window.ElectronCloud.UI = {};

// 初始化UI事件监听
window.ElectronCloud.UI.init = function() {
    const ui = window.ElectronCloud.ui;
    
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
    
    if (ui.sizeSelect) {
        ui.sizeSelect.addEventListener('change', window.ElectronCloud.UI.onSizeChange);
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
    
    // 强制同步一次多选状态
    window.ElectronCloud.UI.onMultiselectToggle();
    
    // 确保初始状态下显示开关的状态正确
    window.ElectronCloud.UI.updateAngular3DToggleState();
    
    // 确保初始状态下清空按钮状态正确
    window.ElectronCloud.UI.updateClearAllSelectionsState();
    
    // 确保初始状态下坐标轴滑动条状态正确
    window.ElectronCloud.UI.updateAxesSizeRangeState();
    
    // 确保初始状态下坐标系正确隐藏（farthestDistance=0时）
    if (window.ElectronCloud.ui.axesSizeRange) {
        window.ElectronCloud.UI.onAxesSizeChange({ target: window.ElectronCloud.ui.axesSizeRange });
    }
    
    console.log('UI 事件监听器初始化完成');
};

// 处理坐标系大小调节
window.ElectronCloud.UI.onAxesSizeChange = function(event) {
    const state = window.ElectronCloud.state;
    const constants = window.ElectronCloud.constants;
    
    if (state.customAxes) {
        const value = parseInt(event.target.value, 10);
        
        // 将滑动条值转换为比例系数 (0-1.05)
        const scaleFactor = value / 100;
        state.axesScaleFactor = scaleFactor;
        
        if (scaleFactor === 0 || state.farthestDistance === 0) {
            // 比例系数为0或轨道半径为0时隐藏坐标系
            state.customAxes.visible = false;
        } else {
            // 显示坐标系并根据比例系数和轨道半径计算大小
            state.customAxes.visible = true;
            
            // 使用实际轨道半径，只有当farthestDistance > 0时才显示坐标系
            const orbitalRadius = Math.max(constants.AXES_BASE_SIZE, state.farthestDistance);
            const targetSize = orbitalRadius * scaleFactor;
            const scale = targetSize / constants.AXES_BASE_SIZE;
            
            state.customAxes.scale.set(scale, scale, scale);
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
    
    if (!ui.sizeSelect) return 0.06;
    
    // 基准点大小（对应基准距离60）
    const v = ui.sizeSelect.value;
    let baseSize;
    if (v === 'small') baseSize = 0.03;
    else if (v === 'medium') baseSize = 0.06;
    else if (v === 'large') baseSize = 0.12;
    else baseSize = 0.06;
    
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
    
    if (state.points && state.points.material && ui.opacityRange) {
        state.points.material.opacity = parseFloat(ui.opacityRange.value);
    }
};

window.ElectronCloud.UI.onCenterLockChange = function() {
    const state = window.ElectronCloud.state;
    const ui = window.ElectronCloud.ui;
    
    if (!state.controls || !ui.centerLock) return;
    
    const locked = ui.centerLock.checked;
    // 禁用平移并锁定目标在原点
    state.controls.enablePan = !locked;
    if (locked) {
        state.controls.target.set(0, 0, 0);
        state.controls.update();
    }
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
    
    // 互斥：取消比照模式
    if (ui.compareToggle && ui.compareToggle.checked) {
        ui.compareToggle.checked = false;
        window.ElectronCloud.UI.onCompareToggle();
    }
    
    const isMultiselect = ui.multiselectToggle && ui.multiselectToggle.checked;
    const label = document.querySelector('label[for="orbital-select"]');
    const controlPanel = document.getElementById('control-panel');
    const multiselectControls = document.querySelector('.multiselect-controls');
    
    if (isMultiselect) {
        if (controlPanel) controlPanel.classList.add('multiselect-active');
        if (multiselectControls) multiselectControls.style.display = 'flex';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = true;
            ui.orbitalSelect.size = 10;
        }
        
        // 多选模式下只禁用3D角向形状开关，保持坐标轴开关可用
        window.ElectronCloud.UI.disableAngular3DToggle();
        
        window.ElectronCloud.UI.updateSelectionCount();
        window.ElectronCloud.UI.refreshSelectStyles();
    } else {
        // 重新启用显示开关功能
        window.ElectronCloud.UI.enableAngular3DToggle();
        
        if (label) label.textContent = '选择轨道';
        if (controlPanel) controlPanel.classList.remove('multiselect-active');
        if (multiselectControls) multiselectControls.style.display = 'none';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = false;
            ui.orbitalSelect.size = 6; // 保持列表样式
        }
        
        // 清除强制样式类
        if (ui.orbitalSelect) {
            const options = ui.orbitalSelect.querySelectorAll('option');
            options.forEach(option => {
                option.classList.remove('force-selected');
            });
        }
        
        // 在退出多选模式时，保留第一个选中的选项，取消其他选项
        if (ui.orbitalSelect && ui.orbitalSelect.selectedOptions.length > 1) {
            const firstSelected = ui.orbitalSelect.selectedOptions[0];
            Array.from(ui.orbitalSelect.options).forEach(option => {
                option.selected = (option === firstSelected);
            });
            const changeEvent = new Event('change', { bubbles: true });
            ui.orbitalSelect.dispatchEvent(changeEvent);
        }
    }
    window.ElectronCloud.UI.updateAngular3DToggleState();
};

// 比照模式功能
window.ElectronCloud.UI.onCompareToggle = function() {
    const ui = window.ElectronCloud.ui;
    
    // 互斥：取消多选模式
    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        ui.multiselectToggle.checked = false;
        window.ElectronCloud.UI.onMultiselectToggle();
    }
    
    const isCompare = ui.compareToggle && ui.compareToggle.checked;
    const label = document.querySelector('label[for="orbital-select"]');
    const controlPanel = document.getElementById('control-panel');
    const multiselectControls = document.querySelector('.multiselect-controls');
    const phaseToggle = document.getElementById('phase-toggle');
    
    if (isCompare) {
        // 禁用渲染相位
        if (phaseToggle) {
            phaseToggle.disabled = true;
            phaseToggle.checked = false;
        }
        
        // 禁用显示开关功能，避免造成bug
        window.ElectronCloud.UI.disableAngular3DToggle();
        
        if (controlPanel) {
            controlPanel.classList.add('multiselect-active');
            controlPanel.classList.add('compare-active');
        }
        if (multiselectControls) multiselectControls.style.display = 'flex';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = true; // 显式开启多选模式
            ui.orbitalSelect.size = 10; // 确保展开显示
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
        
        // 重新启用显示开关功能
        window.ElectronCloud.UI.enableAngular3DToggle();
        
        if (label) label.textContent = '选择轨道';
        if (controlPanel) {
            controlPanel.classList.remove('multiselect-active');
            controlPanel.classList.remove('compare-active');
        }
        if (multiselectControls) multiselectControls.style.display = 'none';
        if (ui.orbitalSelect) {
            ui.orbitalSelect.style.pointerEvents = 'auto';
            ui.orbitalSelect.multiple = false; // 关闭多选模式
            ui.orbitalSelect.size = 6; // 保持列表样式
        }
        
        // 清除强制样式类
        if (ui.orbitalSelect) {
            const options = ui.orbitalSelect.querySelectorAll('option');
            options.forEach(option => {
                option.classList.remove('force-selected');
                option.style.backgroundColor = '';
                option.style.color = '';
            });
        }
        
        // 在退出比照模式时，保留第一个选中的选项，取消其他选项
        if (ui.orbitalSelect && ui.orbitalSelect.selectedOptions.length > 1) {
            const firstSelected = ui.orbitalSelect.selectedOptions[0];
            Array.from(ui.orbitalSelect.options).forEach(option => {
                option.selected = (option === firstSelected);
            });
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
    
    // 多选模式和比照模式的点击处理
    ui.orbitalSelect.addEventListener('mousedown', function(event) {
        const option = event.target;
        if (option.tagName === 'OPTION') {
            // 如果渲染完成且在比照模式下，点击切换轨道可见性
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
                alert('比照模式下最多选择三个轨道');
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
};

// 更新选择计数
window.ElectronCloud.UI.updateSelectionCount = function() {
    const ui = window.ElectronCloud.ui;
    const state = window.ElectronCloud.state;
    
    if (!ui.orbitalSelect) return;
    
    const selectedCount = ui.orbitalSelect.selectedOptions.length;
    const label = document.querySelector('label[for="orbital-select"]');
    
    if (ui.multiselectToggle && ui.multiselectToggle.checked) {
        if (label) label.innerHTML = `选择轨道<br>已选: ${selectedCount}`;
    } else if (ui.compareToggle && ui.compareToggle.checked) {
        // 渲染完成后显示「开关显示」，否则显示「选择轨道」
        if (state.renderingCompleted) {
            if (label) label.innerHTML = `开关显示<br>已选: ${selectedCount}/3`;
        } else {
            if (label) label.innerHTML = `选择轨道<br>已选: ${selectedCount}/3`;
        }
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
        // 清除颜色类，但保留隐藏/可见类（如果渲染完成）
        option.classList.remove('force-selected', 'compare-color-0', 'compare-color-1', 'compare-color-2');
        
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
            option.style.backgroundColor = '';
            option.style.color = '';
        }
    });
};

// 清除所有选择
window.ElectronCloud.UI.clearAllSelections = function() {
    const ui = window.ElectronCloud.ui;
    
    if (!ui.orbitalSelect || !ui.clearAllSelectionsBtn) return;
    
    // 取消所有选中状态
    const options = ui.orbitalSelect.querySelectorAll('option');
    options.forEach(option => {
        option.selected = false;
        option.classList.remove('force-selected');
    });
    
    // 添加清除动画效果
    ui.clearAllSelectionsBtn.style.transform = 'scale(0.95)';
    setTimeout(() => {
        ui.clearAllSelectionsBtn.style.transform = 'scale(1)';
    }, 150);
    
    // 触发change事件更新选择状态
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
    const state = window.ElectronCloud.state;
    const fullscreenBtn = document.getElementById('app-fullscreen-btn');
    
    if (!fullscreenBtn) return;
    
    // 渲染过程中（未完成）禁用全屏按钮
    if (state.isDrawing && !state.samplingCompleted) {
        fullscreenBtn.disabled = true;
        fullscreenBtn.title = '渲染中不可切换全屏';
    } else {
        fullscreenBtn.disabled = false;
        fullscreenBtn.title = document.fullscreenElement ? '退出全屏' : '全屏模式';
    }
};
