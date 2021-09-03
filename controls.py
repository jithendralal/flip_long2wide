import ipywidgets as widgets
from main import IMPLEMENTED


def analysis_type_label():
    text = 'Select Analysis Type:'
    return widgets.HTML(value=f"<h3><font color='teal'>{text}</h3>")


def analysis_types():
    return widgets.ToggleButtons(
        options=IMPLEMENTED,
        description='',
        disabled=False,
        button_style=''
    )


def show_message(text, kind):
    if kind == "success":
        return widgets.HTML(value=f"<h3><font color='green'>{text}</h3>")
    else:
        return widgets.HTML(value=f"<h3><font color='red'>{text}</h3>")
