# VSCode 中 GitHub Copilot Agent 模式免审批功能设置指南

## 概述

GitHub Copilot 的 Agent 模式（也称为 Copilot Edits）允许 AI 直接对代码进行编辑。默认情况下，每次编辑都需要用户手动审批。本指南介绍如何配置免审批（自动批准）功能，以提高开发效率。

---

## 前置要求

1. **VSCode 版本**：需要 VSCode 1.99 或更高版本
2. **GitHub Copilot 扩展**：确保已安装最新版本的 GitHub Copilot 和 GitHub Copilot Chat 扩展
3. **有效订阅**：需要有效的 GitHub Copilot 订阅（Individual、Business 或 Enterprise）

---

## 设置方法

### 方法一：通过 VSCode 设置界面

1. 打开 VSCode
2. 按 `Ctrl + ,`（Windows/Linux）或 `Cmd + ,`（macOS）打开设置
3. 在搜索框中输入 `copilot agent`
4. 找到以下设置项并进行配置：

#### 关键设置项

| 设置项 | 说明 | 推荐值 |
|--------|------|--------|
| `github.copilot.chat.agent.autoApprove` | 自动批准 Agent 的编辑操作 | `true` |
| `github.copilot.chat.agent.runTasks.autoApprove` | 自动批准运行任务 | `true` |
| `github.copilot.chat.agent.runTerminalCommands.autoApprove` | 自动批准终端命令执行 | `false`（安全考虑）|

### 方法二：通过 settings.json 直接配置

1. 按 `Ctrl + Shift + P`（Windows/Linux）或 `Cmd + Shift + P`（macOS）
2. 输入 `Preferences: Open User Settings (JSON)` 并选择
3. 添加以下配置：

```json
{
    "github.copilot.chat.agent.autoApprove": true,
    "github.copilot.chat.agent.runTasks.autoApprove": true,
    "github.copilot.chat.agent.runTerminalCommands.autoApprove": true
}
```

### 方法三：工作区级别设置（推荐用于团队项目）

在项目根目录创建或编辑 `.vscode/settings.json`：

```json
{
    "github.copilot.chat.agent.autoApprove": true,
    "github.copilot.chat.agent.runTasks.autoApprove": true,
    "github.copilot.chat.agent.runTerminalCommands.autoApprove": true
}
```

---

## 详细设置说明

### 1. `github.copilot.chat.agent.autoApprove`

- **作用**：控制是否自动批准 Agent 对文件的编辑操作
- **默认值**：`false`
- **可选值**：
  - `true`：自动批准所有编辑
  - `false`：每次编辑需要手动确认

### 2. `github.copilot.chat.agent.runTasks.autoApprove`

- **作用**：控制是否自动批准运行构建任务（如 `npm build`、`python build_exe.py` 等）
- **默认值**：`false`
- **可选值**：
  - `true`：自动批准任务执行
  - `false`：每次执行需要手动确认

### 3. `github.copilot.chat.agent.runTerminalCommands.autoApprove`

- **作用**：控制是否自动批准在终端执行命令
- **默认值**：`false`
- **可选值**：
  - `true`：自动批准终端命令
  - `false`：每次执行需要手动确认

> ⚠️ **安全提示**：启用终端命令自动批准时请格外小心，因为这可能执行任意命令。建议仅在受信任的项目中启用。

---

## 使用 Agent 模式

### 启动 Agent 模式

1. 按 `Ctrl + Shift + I`（Windows/Linux）或 `Cmd + Shift + I`（macOS）打开 Copilot Chat
2. 或者点击侧边栏的 Copilot 图标
3. 在聊天输入框中，切换到 "Agent" 模式（如果可用）

### 常用 Agent 命令

在 Agent 模式下，可以直接请求代码修改：

```
@workspace 请修复 physics.js 中的波函数计算问题
```

```
@workspace 为 sampling.js 添加单元测试
```

```
@workspace 优化 scene.js 的渲染性能
```

---

## 安全建议

### 推荐的安全配置

对于本项目（氢原子电子云可视化），推荐以下配置：

```json
{
    // 自动批准代码编辑 - 提高开发效率
    "github.copilot.chat.agent.autoApprove": true,
    
    // 自动批准构建任务 - 便于快速测试
    "github.copilot.chat.agent.runTasks.autoApprove": true,
    
    // 终端命令保持手动确认 - 保障安全
    "github.copilot.chat.agent.runTerminalCommands.autoApprove": false
}
```

### 何时应该禁用自动批准

- 处理敏感数据或凭证时
- 在生产环境或生产分支工作时
- 不确定 AI 建议是否正确时
- 进行安全相关的代码修改时

---

## 故障排除

### 问题：设置不生效

1. 确保 VSCode 和 Copilot 扩展已更新到最新版本
2. 重新加载 VSCode 窗口（`Ctrl + Shift + P` → `Developer: Reload Window`）
3. 检查是否有工作区设置覆盖了用户设置

### 问题：Agent 模式不可用

1. 确认 GitHub Copilot Chat 扩展已正确安装
2. 确认已登录 GitHub 账户
3. 确认 Copilot 订阅有效
4. 检查网络连接

### 问题：某些操作仍需要审批

某些高风险操作（如删除文件、执行系统命令）可能仍需要手动确认，这是出于安全考虑的设计。

---

## 相关资源

- [GitHub Copilot 官方文档](https://docs.github.com/en/copilot)
- [VSCode Copilot Chat 文档](https://code.visualstudio.com/docs/copilot/copilot-chat)
- [Copilot 设置参考](https://docs.github.com/en/copilot/configuring-github-copilot/configuring-github-copilot-in-your-environment)

---

## 更新日志

- **2025-12-02**：初始版本，添加 Agent 模式免审批设置指南
