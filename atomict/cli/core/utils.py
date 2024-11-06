def format_pagination_footer(count: int, page_size: int) -> str:
    """Format pagination footer string

    Args:
        count: Total number of items
        page_size: Number of items per page

    Returns:
        Formatted footer string with pagination details
    """
    total_pages = (count - 1) // page_size + 1
    return (
        f"Showing page 1 of {total_pages}, "
        f"Total items: {count}. "
        "Use --all to fetch all results"
    )


def get_pagination_info(results: dict | list) -> tuple[list, str]:
    """Extract pagination information from API response"""
    if isinstance(results, dict):
        items = results.get("results", [])
        count = results.get("count")
        page_size = len(items)

        if count and page_size:
            return items, format_pagination_footer(count, page_size)

        return items, ""

    # All pages were fetched
    return results, f"All {len(results)} results"
