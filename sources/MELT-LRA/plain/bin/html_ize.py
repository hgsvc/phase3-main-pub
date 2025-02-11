#!/usr/bin/env python3

import re
import sys

# ------------------------------------------------------
# main()
# ------------------------------------------------------
def main():

    print("<html>");
    print("<body>");
    print("<pre>");

    for line in sys.stdin:
        line = re.sub(r'<', '&lt;', line.rstrip())
        line = re.sub(r'>', '&gt;', line)
        print(line)
        
    print("</pre>");
    print("</body>");
    print("</html>");

if __name__ == '__main__':
    main()
