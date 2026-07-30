from typing import Optional

def is_not_empty(s: Optional[str]) -> bool:
    return s is not None and s.strip() != ''

def escape_ttl(s: Optional[str]) -> str:
    if s is None:
        return ''
    return (s
            .replace('\\', '\\\\')
            .replace('"', '\\"')
            .replace('\n', '\\n')
            .replace('\r', '\\r')
            .replace('\t', '\\t'))

