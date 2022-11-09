def downLoadjpg(wo='原神芭芭拉', nu=30, html=''):
    from re import findall
    from os import mkdir
    from requests import get
    from getpass import getuser
    from pyautogui import alert
    try:
        mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    except FileExistsError:
        pass
    li = findall(r'https://[^;]+?\.[pngje]{3,4}', html)
    alert('已找到{}张图片'.format(len(li)))
    for i in range(nu):
        try:
            pb = get(li[i]).content
            path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(getuser(), wo, wo, i + 1, nu)
            f = open(path, 'wb')
            f.write(pb)
            f.close()
        except:
            continue
    alert('下载完毕(可能漏下)')

