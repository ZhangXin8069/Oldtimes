def gethtmlbyselenium(
        url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0', ):
    from selenium import webdriver
    from selenium.webdriver.edge.options import Options
    from time import sleep
    # options = Options()
    # options.add_argument('--headless')
    # options.add_argument('--disable-gpu')
    # driver = webdriver.Edge(options=options)
    driver = webdriver.Edge()
    driver.get(url)
    for i in range(5):
        js = "window.scrollTo(0,5000);"
        driver.execute_script(js)
        sleep(0.5)
    html = driver.page_source
    driver.close()
    driver.quit()
    return html


print(gethtmlbyselenium(
    'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=0&latest=0&copyright=0')
)
