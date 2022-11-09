import importlib
import inspect
import os


def filetree(path, depth, f):
    if depth == 0:
        print("文件夹:" + path, file=f)
    for file in os.listdir(path):
        print("|    " * depth + "+--" + file, file=f)
        directory = path + '\\' + file
        if file.endswith(".py"):
            use_inspect(path, depth, f, file)
        if os.path.isdir(directory):
            filetree(directory, depth + 1, f)


def use_inspect(path, depth, f, file):
    fileN = file.replace(".py", "")
    module_name = path.replace(rootdir, "").replace("\\", ".") + "." + fileN
    print("|    " * (depth + 1) + "+--模块" + file, file=f)
    try:
        my_module = importlib.import_module(module_name)
        for name, obj in inspect.getmembers(my_module):
            if inspect.isclass(obj):
                pt0 = "|    " * (depth + 1) + "+--" + str(obj)
                print(pt0, file=f)
                for k, v in obj.__dict__.items():
                    if not (str(k).startswith("_")):
                        pt1 = "|    " * (depth + 2) + "+--" + str(k) + ":" + str(v)
                        print(pt1, file=f)
            if inspect.isfunction(obj):
                pt0 = "|    " * (depth + 1) + "+--" + str(obj)
                print(pt0, file=f)
    except:
        print("导入" + module_name + "失败", file=f)


rootdir = "D:\\Programming\\Python\\Python3.9.10\\Lib\\site-packages\\"
package = "openpyxl"
with open('test.txt', 'w') as f:
    filetree(rootdir + package, 0, f)
