import sys
import os
import socket
import threading
import webbrowser
import mimetypes
from http.server import SimpleHTTPRequestHandler, HTTPServer
import time

# Ensure correct MIME types are known
mimetypes.init()
mimetypes.add_type('application/javascript', '.js')
mimetypes.add_type('application/json', '.json')
mimetypes.add_type('text/css', '.css')
mimetypes.add_type('application/wasm', '.wasm')

def find_free_port():
    """Find an available port on localhost."""
    with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
        s.bind(('', 0))
        _, port = s.getsockname()
    return port

class CORSRequestHandler(SimpleHTTPRequestHandler):
    """Request handler with CORS support and suppressed logging."""
    
    extensions_map = SimpleHTTPRequestHandler.extensions_map.copy()
    extensions_map.update({
        '.js': 'application/javascript',
        '.json': 'application/json',
        '.css': 'text/css',
        '.wasm': 'application/wasm',
    })

    def end_headers(self):
        self.send_header('Access-Control-Allow-Origin', '*')
        super().end_headers()
    
    def log_message(self, format, *args):
        # Suppress logging to keep console clean
        pass

def start_server(port, root_dir):
    """Start the HTTP server serving root_dir on the specified port."""
    # Ensure we are serving files from the correct directory
    os.chdir(root_dir)
    HTTPServer.allow_reuse_address = True
    httpd = HTTPServer(('127.0.0.1', port), CORSRequestHandler)
    httpd.serve_forever()

if __name__ == "__main__":
    # Determine the root directory (handles both script and frozen exe modes)
    if getattr(sys, 'frozen', False):
        # When frozen as exe, use the directory where the exe is located
        # (not sys._MEIPASS which is a temp extraction dir)
        base_path = os.path.dirname(sys.executable)
    else:
        base_path = os.path.dirname(os.path.abspath(__file__))

    # Find a free port
    port = find_free_port()
    
    # Start the server in a daemon thread
    server_thread = threading.Thread(target=start_server, args=(port, base_path))
    server_thread.daemon = True
    server_thread.start()

    # Give server a moment to start
    time.sleep(0.5)

    # Open the default web browser
    url = f"http://127.0.0.1:{port}/index.html"
    webbrowser.open(url)

    # Keep the main thread alive to allow the server thread to run
    try:
        while True:
            time.sleep(1)
    except KeyboardInterrupt:
        sys.exit(0)
