def to_camel(val: str):
    """converts snake_case to camelCase"""
    always_upper = ("sv", "pass", "bed", "re", "ci", "af")
    # ensure all lowercase, split on _
    parts = val.lower().split("_")
    if len(parts) == 1:
        # only one word
        return val

    # first part is lowercase, unless it's an all uppercase word
    camel_parts = [parts[0].upper() if parts[0] in always_upper else parts[0]]
    for snake_part in parts[1:]:
        if snake_part in always_upper:
            # sv_input_file -> SV
            camel_parts.append(snake_part.upper())
        elif camel_parts[-1].isupper():
            # e.g., sv_input_file -> SVinput
            camel_parts.append(snake_part)
        else:
            # sv_input_file -> SVinputFile
            camel_parts.append(snake_part.capitalize())
    return "".join(camel_parts)


def from_camel(val: str):
    """converts camelCase to snake_case"""

    # word boundaries at shift from lower to upper case e.g., camel^Case
    # or upper to lower if several uppercase characters in a row e.g., UPPER^lower
    parts = []
    current_part = None
    for i, c in enumerate(val, 1):
        if current_part is None:
            # beginning of the string
            current_part = c
        else:
            if c.isupper() and current_part[-1].islower():
                # standard break: part^C
                parts.append(current_part)
                current_part = c
            elif len(current_part) == 1:
                # min word length of 2, so always append to current_part if length is 1
                current_part += c
            elif current_part.isupper():
                # current_part length >= 2, all uppercase
                if c.isupper():
                    # more uppercase, append
                    current_part += c
                else:
                    # word break
                    parts.append(current_part)
                    current_part = c
            else:
                # current_part length >= 2, mixed upper/lower or all lower
                if c.isupper():
                    # word break
                    parts.append(current_part)
                    current_part = c
                else:
                    # append and continue
                    current_part += c

        # end of string
        if i == len(val):
            parts.append(current_part)

    return "_".join([x.lower() for x in parts])
