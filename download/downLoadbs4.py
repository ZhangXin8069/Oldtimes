def downLoadbas4(wo='原神芭芭拉', nu=30, html=''):
    from re import findall
    from os import mkdir
    from getpass import getuser
    from pyautogui import alert
    try:
        mkdir(r"C:/Users/{}/Desktop/{}".format(getuser(), wo))
    except FileExistsError:
        pass
    li = findall(r'\"(data:image/.*?)\"', html)
    import base64
    alert('已找到{}张图片'.format(len(li)))
    for i in range(len(li)):
        try:
            da = li[i].split(',')[1]
            pb = base64.b64decode(da)
            path = r'C:/Users/{}/Desktop/{}/{}_{}_{}.png'.format(getuser(), wo, wo, i + 1, nu)
            with open(path, 'wb') as f:
                f.write(pb)
        except Exception as e:
            print(e)
            continue
    alert('下载完毕(可能漏下)')


if __name__ == '__main__':

    def gethtmlbyselenium(
            url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
        from selenium import webdriver
        from time import sleep
        # options = Options()
        # options.add_argument('--headless')
        # options.add_argument('--disable-gpu')
        # driver = webdriver.Edge(options=options)
        driver = webdriver.Edge()
        driver.get(url)
        for i in range(10):
            js = "window.scrollTo(0,10000);"
            driver.execute_script(js)
            print(i)
            sleep(1)
        sleep(60)
        html = driver.page_source
        driver.close()
        driver.quit()
        return html


    url = 'https://lib-lzu.wqxuetang.com/read/pdf?bid=3215540'
    html = gethtmlbyselenium(url)
    downLoadbas4(html=html)
