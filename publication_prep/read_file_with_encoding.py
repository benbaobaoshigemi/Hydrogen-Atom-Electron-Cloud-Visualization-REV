import sys

def read_file(file_path):
    for encoding in ['utf-8', 'gbk']:
        try:
            with open(file_path, 'r', encoding=encoding) as f:
                content = f.read()
                print(f"--- Successfully read with {encoding} ---")
                print(content)
                return
        except UnicodeDecodeError:
            print(f"--- Failed to read with {encoding} ---")
        except Exception as e:
            print(f"--- Error reading with {encoding}: {e} ---")


if __name__ == "__main__":
    if len(sys.argv) > 1:
        read_file(sys.argv[1])
    else:
        print("Please provide a file path.")
