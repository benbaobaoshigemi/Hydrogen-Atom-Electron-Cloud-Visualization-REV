# -*- coding: utf-8 -*-
"""
氢原子电子云可视化 - EXE打包脚本
运行此脚本将项目打包成可独立运行的exe文件
"""

import subprocess
import sys
import os
import shutil

def check_pyinstaller():
    """检查并安装PyInstaller"""
    try:
        import PyInstaller
        print("✓ PyInstaller 已安装")
        return True
    except ImportError:
        print("正在安装 PyInstaller...")
        subprocess.check_call([sys.executable, "-m", "pip", "install", "pyinstaller"])
        print("✓ PyInstaller 安装完成")
        return True

def create_launcher_script():
    """创建启动器脚本"""
    launcher_code = '''# -*- coding: utf-8 -*-
"""
氢原子电子云可视化 - 启动器
双击此exe即可启动程序
"""

import http.server
import socketserver
import webbrowser
import os
import sys
import threading
import time

def get_base_path():
    """获取程序运行的基础路径（资源文件所在目录）"""
    if getattr(sys, 'frozen', False):
        # PyInstaller 打包后，资源文件在 _internal 目录
        exe_dir = os.path.dirname(sys.executable)
        internal_dir = os.path.join(exe_dir, '_internal')
        if os.path.exists(internal_dir) and os.path.exists(os.path.join(internal_dir, 'index.html')):
            return internal_dir
        return exe_dir
    else:
        # 直接运行 Python 脚本时的路径
        return os.path.dirname(os.path.abspath(__file__))

def main():
    # 设置端口
    PORT = 8000
    
    # 切换到资源文件所在目录
    base_path = get_base_path()
    os.chdir(base_path)
    
    print("=" * 50)
    print("   氢原子电子云三维可视化")
    print("   化学教学演示工具")
    print("=" * 50)
    print()
    
    class Handler(http.server.SimpleHTTPRequestHandler):
        def end_headers(self):
            self.send_header('Access-Control-Allow-Origin', '*')
            self.send_header('Cross-Origin-Opener-Policy', 'same-origin')
            self.send_header('Cross-Origin-Embedder-Policy', 'require-corp')
            self.send_header('Cache-Control', 'no-store, no-cache, must-revalidate, max-age=0')
            self.send_header('Pragma', 'no-cache')
            self.send_header('Expires', '0')
            super().end_headers()
        
        def log_message(self, format, *args):
            # 简化日志输出
            pass
    
    # 寻找可用端口
    while True:
        try:
            httpd = socketserver.TCPServer(("", PORT), Handler)
            break
        except OSError:
            PORT += 1
            if PORT > 8100:
                print("错误：无法找到可用端口")
                input("按回车键退出...")
                sys.exit(1)
    
    print(f"✓ 服务已启动")
    print(f"✓ 访问地址: http://localhost:{PORT}")
    print()
    print("提示：")
    print("  - 请勿关闭此窗口")
    print("  - 浏览器会自动打开")
    print("  - 关闭此窗口即可退出程序")
    print()
    print("-" * 50)
    
    # 延迟打开浏览器，确保服务器已完全启动
    def open_browser():
        time.sleep(1)
        webbrowser.open(f"http://localhost:{PORT}")
    
    browser_thread = threading.Thread(target=open_browser)
    browser_thread.daemon = True
    browser_thread.start()
    
    try:
        httpd.serve_forever()
    except KeyboardInterrupt:
        print()
        print("程序已退出")
    finally:
        httpd.server_close()

if __name__ == "__main__":
    main()
'''
    
    with open("launcher.py", "w", encoding="utf-8") as f:
        f.write(launcher_code)
    print("✓ 启动器脚本已创建")

def build_exe():
    """使用PyInstaller打包"""
    print()
    print("开始打包...")
    print()
    
    # 获取当前目录
    current_dir = os.path.dirname(os.path.abspath(__file__))
    
    # 需要打包的数据文件
    data_files = [
        "index.html",
        "style.css",
        "main.js",
        "core.js",
        "orbital.js",
        "physics.js",
        "sampling.js",
        "sampling-worker.js",
        "scene.js",
        "ui.js",
        "visualization.js",
        "data_panel.js",
        "gesture_v2.js",
        "gesture_config.json",
    ]
    
    # 构建 --add-data 参数
    add_data_args = []
    for f in data_files:
        if os.path.exists(f):
            add_data_args.extend(["--add-data", f"{f};."])
    
    # PyInstaller 命令
    cmd = [
        sys.executable, "-m", "PyInstaller",
        "--onedir",  # 打包成一个文件夹
        "--windowed",  # 不显示控制台窗口（但我们需要显示，所以改用console）
        "--console",  # 显示控制台窗口
        "--name", "氢原子电子云可视化",
        "--noconfirm",  # 不询问确认
        "--clean",  # 清理临时文件
    ]
    
    # 添加数据文件
    cmd.extend(add_data_args)
    
    # 添加启动器脚本
    cmd.append("launcher.py")
    
    print("执行命令:")
    print(" ".join(cmd[:10]) + " ...")
    print()
    
    # 执行打包
    result = subprocess.run(cmd, capture_output=False)
    
    if result.returncode == 0:
        print()
        print("=" * 50)
        print("✓ 打包成功！")
        print("=" * 50)
        print()
        print(f"输出目录: {os.path.join(current_dir, 'dist', '氢原子电子云可视化')}")
        print()
        print("使用说明:")
        print("1. 将 dist/氢原子电子云可视化 整个文件夹复制到目标电脑")
        print("2. 双击文件夹中的 '氢原子电子云可视化.exe' 即可运行")
        print()
    else:
        print()
        print("✗ 打包失败，请检查错误信息")
        
    return result.returncode == 0

def cleanup():
    """清理临时文件"""
    temp_files = ["launcher.py", "launcher.spec"]
    for f in temp_files:
        if os.path.exists(f):
            try:
                os.remove(f)
            except:
                pass
    
    # 清理 build 目录
    if os.path.exists("build"):
        try:
            shutil.rmtree("build")
        except:
            pass

def main():
    print()
    print("=" * 50)
    print("   氢原子电子云可视化 - EXE打包工具")
    print("=" * 50)
    print()
    
    # 切换到脚本所在目录
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    
    # 检查并安装 PyInstaller
    check_pyinstaller()
    
    # 创建启动器脚本
    create_launcher_script()
    
    # 打包
    success = build_exe()
    
    # 清理临时文件
    cleanup()
    
    print()
    if success:
        input("按回车键退出...")
    else:
        input("打包失败，按回车键退出...")

if __name__ == "__main__":
    main()
