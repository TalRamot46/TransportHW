"""Logged tables and headings, so the reports carry no format strings of their own."""

import logging

logger = logging.getLogger(__name__)

def log_section(title, *notes):
    """Logs a question heading followed by any explanatory lines."""
    logger.info(f"\n=== {title} ===")
    for note in notes:
        logger.info(note)

def log_table(headers, rows, fmt='{:.6f}'):
    """Logs a table; row values are formatted with `fmt` unless already strings."""
    cells = [[v if isinstance(v, str) else fmt.format(v) for v in row] for row in rows]
    widths = [max([len(h)] + [len(row[i]) for row in cells])
              for i, h in enumerate(headers)]
    rule = '-' * (sum(widths) + 3 * len(widths))

    logger.info(rule)
    logger.info(' | '.join(h.ljust(w) for h, w in zip(headers, widths)))
    logger.info(rule)
    for row in cells:
        logger.info(' | '.join(v.ljust(w) for v, w in zip(row, widths)))
    logger.info(rule)
