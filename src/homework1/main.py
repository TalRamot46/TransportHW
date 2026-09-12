"""Entry point: runs every question's checks and writes its figures."""

import logging
from homework1 import q1, q2, q3, q4, q5
from homework1.figures import figs_dir

logging.basicConfig(level=logging.INFO, format='%(message)s')
logger = logging.getLogger(__name__)

def main():
    """Runs Questions 1 to 5 in order."""
    figs = figs_dir()
    for question in (q1, q2, q3, q4, q5):
        question.report(figs)
    logger.info(f"\nAll figures written to {figs}")

if __name__ == "__main__":
    main()
