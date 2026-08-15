"""Entry point: runs Questions 1, 3, 4 and 5 and writes their figures."""

import logging
from homework3 import q1, q3, q4, q5
from homework3.figures import figs_dir

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def main():
    """Runs the solved questions in order."""
    figs = figs_dir()
    for question in (q1, q3, q4, q5):
        question.report(figs)
    logger.info(f"\nAll figures written to {figs}")

if __name__ == "__main__":
    main()
