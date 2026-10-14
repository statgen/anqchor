"""Allow ``python -m anqchor``; the generated Makefiles use this form so that step
commands run with the same interpreter and environment as the parent invocation."""

from anqchor.cli.main import main

if __name__ == "__main__":
    main()
