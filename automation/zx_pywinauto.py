from pywinauto.application import Application

# Popen('calc.exe', shell=True)
app = Application()
app.start(("notepad.exe"))
# app = Application().connect(path=r"c:\windows\system32\notepad.exe")
# dlg_spec = app.UntitledNotepad
# app.UntitledNotepad.type_keys("%FX")
# app.UntitledNotepad.menu_select("File->SaveAs")
# app.SaveAs.ComboBox5.select("UTF-8")
# app.SaveAs.edit1.set_text("Example-utf8.txt")
# app.SaveAs.Save.click()
# dlg_spec = app.window(title_re='.* - Notepad$').window(class_name='Edit')
dlg_spec = app.window(title='无标题 - 记事本')
dlg_spec.minimize()
