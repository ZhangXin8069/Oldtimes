

def input_set(a={'a': 5, 'b': '6', 'c': 9}):
    from pyautogui import prompt, alert
    for ai in a:
        ay = a[ai]
        t = type(a[ai])
        while True:
            a[ai] = ay
            try:
                p = prompt("{}|默认:{};类型:{};修改:输入'@'开始修改;退出:输入'#'准备退出".format(ai, a[ai], t))
                if p == '':
                    break
                if p == '#':
                    exit(0)
                if p == '@':
                    a[ai] = prompt(ai)
                    if type(eval(a[ai])) != t:
                        alert('输入有误')
                    else:
                        alert('输入完成')
                        break
            except Exception as e:
                alert(e)
                alert('重新输入')
    return a




