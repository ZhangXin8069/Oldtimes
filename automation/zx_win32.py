import time

import win32api
import win32con
import win32gui


def open_max_name(path="notepad.exe", name=u'无标题 - 记事本'):
    win32api.ShellExecute(0, 'open', path, '', '', 1)
    time.sleep(0.5)
    hwnd = win32gui.FindWindow(None, name)
    print(hwnd, type(hwnd))
    win32gui.SetForegroundWindow(hwnd)
    win32gui.PostMessage(hwnd, win32con.WM_SYSCOMMAND, win32con.SC_MAXIMIZE, 0)


def open_max_hwnd(hwnd=65539):
    print(hwnd, type(hwnd))
    win32gui.SetForegroundWindow(hwnd)
    win32gui.PostMessage(hwnd, win32con.WM_SYSCOMMAND, win32con.SC_MAXIMIZE, 0)


def get_all_windows():
    hWnd_list = []
    win32gui.EnumWindows(lambda hWnd, param: param.append(hWnd), hWnd_list)
    # print(hWnd_list)
    return hWnd_list


def get_son_windows(parent):
    hWnd_child_list = []
    win32gui.EnumChildWindows(parent, lambda hWnd, param: param.append(hWnd), hWnd_child_list)
    print(hWnd_child_list)
    return hWnd_child_list


def get_title(hwnd):
    title = win32gui.GetWindowText(hwnd)
    print('窗口标题:%s' % (title))
    return title


def get_clasname(hwnd):
    clasname = win32gui.GetClassName(hwnd)
    print('窗口类名:%s' % (clasname))
    return clasname


def set_top(hwnd):
    win32gui.SetWindowPos(hwnd, win32con.HWND_TOPMOST, 0, 0, 0, 0,
                          win32con.SWP_NOMOVE | win32con.SWP_NOACTIVATE | win32con.SWP_NOOWNERZORDER | win32con.SWP_SHOWWINDOW | win32con.SWP_NOSIZE)


def set_down(hwnd):
    win32gui.SetWindowPos(hwnd, win32con.HWND_NOTOPMOST, 0, 0, 0, 0,
                          win32con.SWP_SHOWWINDOW | win32con.SWP_NOSIZE | win32con.SWP_NOMOVE)


# 根据窗口名称获取句柄
def get_hwnd_from_name(name):
    hWnd_list = get_all_windows()
    for hwd in hWnd_list:
        title = get_title(hwd)
        if title == name:
            return hwd


def xianshi(name):
    hwd = get_hwnd_from_name(name)
    win32gui.ShowWindow(hwd, win32con.SW_SHOW)


def yingcang(name):
    hwd = get_hwnd_from_name(name)
    win32gui.ShowWindow(hwd, win32con.SW_HIDE)


# 获取右下角托盘的任务句柄
def get_tuopan_hwd():
    handle = win32gui.FindWindow("Shell_TrayWnd", None)
    hWnd_child_list = get_son_windows(handle)[1:]
    tuopan_hwd_list = []
    flag = False
    for i in hWnd_child_list:
        if get_clasname(i) == 'TrayNotifyWnd':
            flag = True
        if flag:
            tuopan_hwd_list.append(i)
    return tuopan_hwd_list


#
# if __name__ == '__main__':
#     # dic = {win32gui.GetWindowText(i): i for i in get_all_windows()}
#     # # dic = {dic[i]: i for i in dic.keys()}
#     # print(dic)
#
if __name__ == '__main__':
    # open_max_name()
    # open_max_hwnd(920174)
    dic = {win32gui.GetWindowText(i): i for i in get_all_windows()}
    print(dic)
