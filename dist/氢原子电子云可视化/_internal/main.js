// 主入口文件 - 启动应用程序
// 依赖模块加载顺序：
// 1. core.js - 核心变量和状态管理
// 2. scene.js - Three.js 场景管理
// 3. orbital.js - 轨道管理和数据处理
// 4. sampling.js - 点采样逻辑
// 5. visualization.js - 3D 可视化效果
// 6. ui.js - UI 控件交互
// 7. main.js - 主入口（本文件）

// 确保所有依赖模块都已加载
function checkDependencies() {
    const requiredModules = [
        'ElectronCloud',
        'ElectronCloud.Scene',
        'ElectronCloud.Orbital',
        'ElectronCloud.Sampling',
        'ElectronCloud.Visualization',
        'ElectronCloud.UI'
    ];
    
    for (const module of requiredModules) {
        const parts = module.split('.');
        let obj = window;
        for (const part of parts) {
            if (!obj[part]) {
                console.error(`依赖模块缺失: ${module}`);
                return false;
            }
            obj = obj[part];
        }
    }
    
    return true;
}

// 应用程序启动函数
function startApplication() {
    console.log('开始启动 ElectronCloud 应用程序...');
    
    // 检查依赖
    if (!checkDependencies()) {
        console.error('依赖检查失败，无法启动应用程序');
        return;
    }
    
    // 检查必要的外部库
    if (typeof THREE === 'undefined') {
        console.error('Three.js 库未加载');
        return;
    }
    
    if (typeof window.Hydrogen === 'undefined') {
        console.error('Hydrogen 物理模型库未加载');
        return;
    }
    
    try {
        // 初始化应用程序
        window.ElectronCloud.init();
        console.log('ElectronCloud 应用程序启动成功！');
    } catch (error) {
        console.error('应用程序启动失败:', error);
    }
}

// 等待页面加载完成后启动
if (document.readyState === 'loading') {
    document.addEventListener('DOMContentLoaded', startApplication);
} else {
    // 如果页面已经加载完成，立即启动
    startApplication();
}

// 导出启动函数供外部调用（如果需要）
window.startElectronCloud = startApplication;
