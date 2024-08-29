#!/usr/bin/env python
"""Django's command-line utility for administrative tasks."""
import os
import sys


import os
import sys

def main():
    """Run administrative tasks."""
    os.environ.setdefault('DJANGO_SETTINGS_MODULE', 'compound_management.settings')
    try:
        # Ensure that Django is installed and can be imported
        import django
        from django.core.management import execute_from_command_line
    except ImportError as exc:
        # Provide a more detailed error message to help with troubleshooting
        raise ImportError(
            "Couldn't import Django. Make sure it's installed and "
            "available on your PYTHONPATH environment variable. "
            "Ensure your virtual environment is activated."
        ) from exc
    except Exception as e:
        # Catch other potential exceptions and provide a meaningful message
        print(f"An unexpected error occurred: {e}")
        sys.exit(1)

    execute_from_command_line(sys.argv)


if __name__ == '__main__':
    main()
