from rich.table import Table
from rich.box import ASCII_DOUBLE_HEAD
from typing import List, Tuple, Any, Optional, Callable


table_0 = Table(
    title_style="bold",
    title_justify="center",
    box=ASCII_DOUBLE_HEAD,
    caption_justify="center",
    show_header=True,
    header_style="bold cyan",
    # show_footer=True,
    # footer_style="bold cyan",
    highlight=True,
)


def create_table(
    columns: List[Tuple[str, str, Optional[Callable]]], 
    items: List[dict],
    *,
    title: Optional[str] = None,
    box: Any = ASCII_DOUBLE_HEAD,
    header_style: str = "bold cyan",
    title_style: str = "bold",
    caption: Optional[str] = None,
) -> Table:
    """Create a styled Rich table with the given columns and items.
    
    Args:
        columns: List of (header, key, formatter) tuples. 
                formatter is optional and defaults to str
        items: List of dictionaries containing the data
        title: Optional table title
        box: Rich box style to use
        header_style: Style for header row
        title_style: Style for title
        footer: Optional footer text
    """
    table = Table(
        title=title,
        title_style=title_style,
        title_justify="center",
        box=box,
        caption_justify="center",
        show_header=True,
        header_style=header_style,
        highlight=True,
    )

    # Add columns
    for header, _, _ in columns:
        table.add_column(header)

    # Add rows
    for item in items:
        row = []
        for header, key, formatter in columns:
            value = item.get(key, "")
            # Use str as default formatter if none provided
            format_func = formatter if formatter is not None else str
            row.append(format_func(value) if value is not None else "")
        table.add_row(*row)

    if caption:
        table.caption = caption

    return table
