import os
import glob
import sys
import re

if len(sys.argv) != 3:
    print("Usage: rename_by_prefix OLD_CONTENT NEW_CONTENT")
else:
    old_content = sys.argv[1]
    new_content = sys.argv[2]

    for filename in glob.glob("*" + old_content + "*"):
        os.rename(filename, re.sub(old_content, new_content, filename))
