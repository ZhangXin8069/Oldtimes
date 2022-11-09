import pyautogui
# # pyautogui.FAILSAFE = True
# # pyautogui.PAUSE = 1
# # print(pyautogui.size())   # 返回所用显示器的分辨率； 输出：Size(width=1920, height=1080)
# # width, height = pyautogui.size()
# # print(width, height)  # 1920 1080
# # # pyautogui.moveTo(100,300,duration=1)
# # pyautogui.moveRel(100,500,duration=4)   # 第一个参数是左右移动像素值，第二个是上下，
# # print(pyautogui.position())   # 得到当前鼠标位置；输出：Point(x=200, y=800)
# # # # 点击鼠标
# # # pyautogui.click(10,10)   # 鼠标点击指定位置，默认左键
# # # pyautogui.click(10,10,button='left')  # 单击左键
# # # pyautogui.click(1000,300,button='right')  # 单击右键
# # # pyautogui.click(1000,300,button='middle')  # 单击中间
# pyautogui.doubleClick(10,10)  # 指定位置，双击左键
# pyautogui.rightClick(10,10)   # 指定位置，双击右键
# pyautogui.middleClick(10,10)  # 指定位置，双击中键
# pyautogui.mouseDown()   # 鼠标按下
# pyautogui.mouseUp()    # 鼠标释放
#
# pyautogui.dragRel(100,500,duration=4)   # 第一个参数是左右移动像素值，第二个是上下，