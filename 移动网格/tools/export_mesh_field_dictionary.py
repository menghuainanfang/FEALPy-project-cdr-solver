"""Write the human-readable field table from the UI's shared dictionary."""

from pathlib import Path
from src.mesh_trace import FIELD_DICTIONARY


def main():
    data = FIELD_DICTIONARY
    text = f"# 网格诊断字段词典\n\n版本：`{data['version']}`。定义源为 `src/mesh_field_dictionary.json`；E/H面板问号使用同一份定义。\n\n{data['index_policy']}\n\n{data['alias_scope']}\n\n"
    text += "|字段|UI名称|单位|定义|历史别名|\n|---|---|---|---|---|\n"
    for field in data["fields"]:
        text += "|"+"|".join((field["key"], field["label"], field["unit"],
                               field["definition"].replace("|", "\\|"), ", ".join(field["aliases"]) or "无"))+"|\n"
    path = Path(__file__).resolve().parents[1] / "docs/mentor_2026-10-23/mesh_field_dictionary_v2.md"
    path.write_text(text, encoding="utf-8")
    print(path)


if __name__ == "__main__":
    main()
