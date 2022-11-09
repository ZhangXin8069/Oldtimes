from __future__ import annotations

import re

from manimlib.mobject.svg.string_mobject import StringMobject
from manimlib.utils.tex_file_writing import display_during_execution
from manimlib.utils.tex_file_writing import tex_content_to_svg_file

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from colour import Color
    import re
    from typing import Iterable, Union

    from manimlib.mobject.types.vectorized_mobject import VGroup

    ManimColor = Union[str, Color]
    Span = tuple[int, int]
    Selector = Union[
        str,
        re.Pattern,
        tuple[Union[int, None], Union[int, None]],
        Iterable[Union[
            str,
            re.Pattern,
            tuple[Union[int, None], Union[int, None]]
        ]]
    ]


SCALE_FACTOR_PER_FONT_POINT = 0.001


class MTex(StringMobject):
    CONFIG = {
        "font_size": 48,
        "alignment": "\\centering",
        "tex_environment": "align*",
        "tex_to_color_map": {},
        "template": "",
        "additional_preamble": "",
    }

    def __init__(self, tex_string: str, **kwargs):
        # Prevent from passing an empty string.
        if not tex_string.strip():
            tex_string = "\\\\"
        self.tex_string = tex_string
        super().__init__(tex_string, **kwargs)

        self.set_color_by_tex_to_color_map(self.tex_to_color_map)
        self.scale(SCALE_FACTOR_PER_FONT_POINT * self.font_size)

    @property
    def hash_seed(self) -> tuple:
        return (
            self.__class__.__name__,
            self.svg_default,
            self.path_string_config,
            self.base_color,
            self.isolate,
            self.protect,
            self.tex_string,
            self.alignment,
            self.tex_environment,
            self.tex_to_color_map,
            self.template,
            self.additional_preamble
        )

    def get_file_path_by_content(self, content: str) -> str:
        with display_during_execution(f"Writing \"{self.tex_string}\""):
            file_path = tex_content_to_svg_file(
                content, self.template, self.additional_preamble
            )
        return file_path

    # Parsing

    @staticmethod
    def get_command_matches(string: str) -> list[re.Match]:
        # Lump together adjacent brace pairs
        pattern = re.compile(r"""
            (?P<command>\\(?:[a-zA-Z]+|.))
            |(?P<open>{+)
            |(?P<close>}+)
        """, flags=re.X | re.S)
        result = []
        open_stack = []
        for match_obj in pattern.finditer(string):
            if match_obj.group("open"):
                open_stack.append((match_obj.span(), len(result)))
            elif match_obj.group("close"):
                close_start, close_end = match_obj.span()
                while True:
                    if not open_stack:
                        raise ValueError("Missing '{' inserted")
                    (open_start, open_end), index = open_stack.pop()
                    n = min(open_end - open_start, close_end - close_start)
                    result.insert(index, pattern.fullmatch(
                        string, pos=open_end - n, endpos=open_end
                    ))
                    result.append(pattern.fullmatch(
                        string, pos=close_start, endpos=close_start + n
                    ))
                    close_start += n
                    if close_start < close_end:
                        continue
                    open_end -= n
                    if open_start < open_end:
                        open_stack.append(((open_start, open_end), index))
                    break
            else:
                result.append(match_obj)
        if open_stack:
            raise ValueError("Missing '}' inserted")
        return result

    @staticmethod
    def get_command_flag(match_obj: re.Match) -> int:
        if match_obj.group("open"):
            return 1
        if match_obj.group("close"):
            return -1
        return 0

    @staticmethod
    def replace_for_content(match_obj: re.Match) -> str:
        return match_obj.group()

    @staticmethod
    def replace_for_matching(match_obj: re.Match) -> str:
        if match_obj.group("command"):
            return match_obj.group()
        return ""

    @staticmethod
    def get_attr_dict_from_command_pair(
        open_command: re.Match, close_command: re.Match
    ) -> dict[str, str] | None:
        if len(open_command.group()) >= 2:
            return {}
        return None

    def get_configured_items(self) -> list[tuple[Span, dict[str, str]]]:
        return [
            (span, {})
            for selector in self.tex_to_color_map
            for span in self.find_spans_by_selector(selector)
        ]

    @staticmethod
    def get_color_command(rgb_hex: str) -> str:
        rgb = MTex.hex_to_int(rgb_hex)
        rg, b = divmod(rgb, 256)
        r, g = divmod(rg, 256)
        return f"\\color[RGB]{{{r}, {g}, {b}}}"

    @staticmethod
    def get_command_string(
        attr_dict: dict[str, str], is_end: bool, label_hex: str | None
    ) -> str:
        if label_hex is None:
            return ""
        if is_end:
            return "}}"
        return "{{" + MTex.get_color_command(label_hex)

    def get_content_prefix_and_suffix(
        self, is_labelled: bool
    ) -> tuple[str, str]:
        prefix_lines = []
        suffix_lines = []
        if not is_labelled:
            prefix_lines.append(self.get_color_command(
                self.color_to_hex(self.base_color)
            ))
        if self.alignment:
            prefix_lines.append(self.alignment)
        if self.tex_environment:
            prefix_lines.append(f"\\begin{{{self.tex_environment}}}")
            suffix_lines.append(f"\\end{{{self.tex_environment}}}")
        return (
            "".join([line + "\n" for line in prefix_lines]),
            "".join(["\n" + line for line in suffix_lines])
        )

    # Method alias

    def get_parts_by_tex(self, selector: Selector) -> VGroup:
        return self.select_parts(selector)

    def get_part_by_tex(self, selector: Selector, **kwargs) -> VGroup:
        return self.select_part(selector, **kwargs)

    def set_color_by_tex(self, selector: Selector, color: ManimColor):
        return self.set_parts_color(selector, color)

    def set_color_by_tex_to_color_map(
        self, color_map: dict[Selector, ManimColor]
    ):
        return self.set_parts_color_by_dict(color_map)

    def get_tex(self) -> str:
        return self.get_string()


class MTexText(MTex):
    CONFIG = {
        "tex_environment": None,
    }
