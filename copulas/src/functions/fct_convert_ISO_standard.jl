function convert_ISO_standard(string)
    """
    Converts string to ISO standard. Timezone part dropped

    Input: 
        - string: string formatted as y-m-d HH:MM:SSzzzz
    Output:
        - DateTime object in the same length as string
    """
    date_format = DateFormat("y-m-d HH:MM:SSzzzz")
    x = DateTime.(string, date_format)
    return (x)
end