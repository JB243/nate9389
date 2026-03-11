import io
import requests
import pandas as pd
import numpy as np
import time
from bs4 import BeautifulSoup, Tag, NavigableString
import html2text
import re
import argparse

# Ensure there's no more than one consecutive newline
def format_output(outputs):
    # Join the outputs and ensure there's exactly one newline between non-empty lines
    return '\n\n'.join(line for line in outputs if line)

def html_to_markdown(html_content):
    # Create an html2text converter object
    h = html2text.HTML2Text()
    h.ignore_links = False

    # Remove occurrences of the strings "[본문]" and "[목차]"
    def remove_unwanted_strings(text):
        text = text.replace('[본문]', '')
        text = text.replace('[목차]', '')
        return text

    # Convert HTML to markdown
    def convert_to_markdown(item):
        if isinstance(item, (Tag, NavigableString)):
            item = str(item)
        markdown = h.handle(item)
        
        # Convert the domain as specified
        markdown = markdown.replace('https://nate9389.tistory.com/', 'https://jb243.github.io/pages/')
        markdown = markdown.replace('1st', '1<sup>st</sup>')
        markdown = markdown.replace('2nd', '2<sup>nd</sup>')
        markdown = markdown.replace('3rd', '3<sup>rd</sup>')
        markdown = markdown.replace('4th', '4<sup>th</sup>')
        markdown = markdown.replace('5th', '5<sup>th</sup>')
        markdown = markdown.replace('6th', '6<sup>th</sup>')
        markdown = markdown.replace('7th', '7<sup>th</sup>')
        markdown = markdown.replace('8th', '8<sup>th</sup>')
        markdown = markdown.replace('9th', '9<sup>th</sup>')
        markdown = markdown.replace('10th', '10<sup>th</sup>')
        markdown = markdown.replace('11th', '11<sup>th</sup>')
        markdown = markdown.replace('12th', '12<sup>th</sup>')
        markdown = markdown.replace('13th', '13<sup>th</sup>')
        markdown = remove_unwanted_strings(markdown).strip()
        
        return markdown

    # If it's a list (ResultSet), iterate and convert each item
    if isinstance(html_content, list):
        markdown_contents = []
        for item in html_content:
            converted_md = convert_to_markdown(item)
            if converted_md:
                prefix = detect_margin_and_prepend(item)
                markdown_contents.append(prefix + converted_md)
        return markdown_contents

    # For single items
    if isinstance(html_content, (Tag, NavigableString)):
        prefix = detect_margin_and_prepend(html_content)
        return prefix + convert_to_markdown(html_content)

    return ""

def detect_margin_and_prepend(parsed_item):
    # If the parsed item is effectively empty, return an empty string
    if not parsed_item.get_text(strip=True):
        return ''

    style = parsed_item.get('style', '')

    # Detect the number of 'em' in margin-left using regex
    margin_left_matches = re.findall(r'margin-left:\s*(\d+)em', style)
    shorthand_margin_matches = re.findall(r'margin:\s*\d+px\s*\d+px\s*\d+px\s*(\d+)em', style)

    # Aggregate both types of matches
    all_matches = margin_left_matches + shorthand_margin_matches

    # If there's a match and it's a valid number, calculate the prefix
    if all_matches and all_matches[0].isdigit():
        em_increment = 2  # Each 'em' increment adds a '>'
        number_of_em = int(all_matches[0])
        prefix_count = number_of_em // em_increment
        prefix = '>' * (prefix_count-1) + ' '
    else:
        prefix = ''

    return prefix

def add_blank_lines(markdown_list):
    return [item + '\n' for item in markdown_list]


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--idx', type=int, help='Input idx')
    args = parser.parse_args()

    url = 'https://nate9389.tistory.com/' + str(args.idx)
    github_session = requests.Session()
    download = github_session.get(url).content
    html = download
    
    # Parse the HTML with BeautifulSoup
    soup = BeautifulSoup(html, 'html.parser')

    # Find the desired content by its attributes
    size_pattern = re.compile(r'size\d+')
    parsed_items = soup.find_all('p', {'data-ke-size': size_pattern})
    
    markdown_outputs = html_to_markdown(parsed_items)
    markdown_outputs[0] = '## ' + markdown_outputs[0]
    
    markdown_outputs = add_blank_lines(markdown_outputs)

    for output in markdown_outputs:
        print(output)
