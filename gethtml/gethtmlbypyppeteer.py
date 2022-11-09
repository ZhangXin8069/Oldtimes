def gethtmlbypyppteer(
        url=r'https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=1&latest=0&copyright=0',):
    import asyncio
    from pyppeteer import launch
    from time import sleep
    async def asgethtml(url):
        browser = await launch({
            'headless': False,
            'devtools': True,
            'args': [
                '--disable-extensions',
                '--hide-scrollbars',
                '--disable-bundled-ppapi-flash',
                '--mute-audio',
                '--no-sandbox',
                '--disable-setuid-sandbox',
                '--disable-gpu',
            ],
            'dumpio': True,
        })
        page = await browser.newPage()
        await page.setUserAgent(
            "Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/97.0.4692.99 Safari/537.36 Edg/97.0.1072.69"
        )
        await page.evaluateOnNewDocument('()=>{Object.defineProperties(navigator,{webdriver:{get:()=>false}})}')
        await page.goto(url)
        for i in range(5):
            await page.evaluate('window.scrollTo(0,5000)')
            sleep(0.5)
        text = await page.content()
        await browser.close()
        return text

    m = asyncio.ensure_future(asgethtml(url))
    asyncio.get_event_loop().run_until_complete(m)
    return m.result()
def F(x):
    return x**x

__all__='gethtmlbypyppteer'
if __name__=='__main__':
    print(gethtmlbypyppteer('https://image.baidu.com/search/index?ct=201326592&z=9&tn=baiduimage&ipn=r&word=%E5%8E%9F%E7%A5%9E%E8%8A%AD%E8%8A%AD%E6%8B%89&pn=0&istype=2&ie=utf-8&oe=utf-8&cl=2&lm=-1&st=-1&fr=&fmq=&ic=0&se=&sme=&width=0&height=0&face=0&hd=0&latest=0&copyright=0'))
    # print(gethtmlbypyppteer())