def downLoadwebp(wo='原神芭芭拉', nu=30, html=''):
    from re import findall
    from os import mkdir
    from io import BytesIO
    from PIL import Image
    from requests import get
    from getpass import getuser
    from pyautogui import alert
    try:
        mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    except FileExistsError:
        pass
    li = findall(r'\"(https://[^;]+?)\"', html)
    li = list(set(li))
    alert('已找到{}张图片'.format(len(li)))
    for i in range(nu):
        try:
            pb = BytesIO(get(li[i]).content)
            print(li[i])
            im = Image.open(pb)
            if im.mode == "RGBA":
                im.load()  # required for png.split()
                background = Image.new("RGB", im.size, (255, 255, 255))
                background.paste(im, mask=im.split()[3])
            path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(getuser(), wo, wo, i + 1, nu)
            im.save(path, 'JPEG')
        except Exception as e:
            print(e)
            continue
    alert('下载完毕(可能漏下)')

