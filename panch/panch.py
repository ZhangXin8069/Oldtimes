from time import sleep

from pyautogui import hotkey, keyDown, keyUp, locateCenterOnScreen, click, locateOnScreen, doubleClick, moveTo
from pyperclip import copy


def uwxn():
    sleep(2)
    keyDown('ctrl')
    hotkey('r')
    keyUp('ctrl')


def enterUrl(url=r'http://my.lzu.edu.cn:8080/login?service=http://my.lzu.edu.cn'):
    keyDown('win')
    hotkey('r')
    keyUp('win')
    copy(url)
    keyDown('ctrl')
    hotkey('v')
    keyUp('ctrl')

    hotkey('enter')


def clickimg(img=r"..\img\max.png"):
    moveTo(10, 10)
    sleep(2)
    click(locateCenterOnScreen(img, confidence=0.8, grayscale=True, region=(0, 0, 1920, 1080)))


def enterHtml():
    from pyautogui import hotkey
    uwxn()
    for i in range(5):
        hotkey('tab')
    hotkey('enter')
    uwxn()
    for i in range(5):
        hotkey('tab')
    hotkey('down')
    hotkey('enter')
    hotkey('enter')


def enterWeb():
    keyDown('alt')
    hotkey('tab')
    clickimg(r"..\img\ljvbdaxt.png")

    keyUp('alt')


def max():
    try:
        t = locateOnScreen(r"..\img\max1.png")
        t = (t[0] + t[3], t[1])
        doubleClick(t)
        t = locateOnScreen(r"..\img\max2.png")
        t = (t[0] + t[3], t[1])
        doubleClick(t)
    except TypeError:
        pass
    clickimg(r"..\img\max3.png")
    clickimg(r"..\img\max4.png")


def clicks():
    clickimg(r"..\img\jxkh.png")
    clickimg(r"..\img\uhbk.png")
    clickimg(r"..\img\qtd;.png")


def panch():
    enterUrl()
    enterWeb()
    max()
    enterHtml()
    clicks()


if __name__ == '__main__':
    # enterUrl()
    # enterHtml()
    panch()
