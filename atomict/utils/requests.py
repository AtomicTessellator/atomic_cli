import requests
from tqdm import tqdm


def download_large_file(url, destination_file):
    # Send an HTTP request to get the file size
    response = requests.get(url, stream=True)
    total_size = int(response.headers.get('content-length', 0))

    # Initialize tqdm progress bar with the total file size
    progress_bar = tqdm(total=total_size, unit='iB', unit_scale=True)

    with open(destination_file, 'wb') as file:
        for chunk in response.iter_content(chunk_size=8192):
            # Write the chunk to the file and update the progress bar
            file.write(chunk)
            progress_bar.update(len(chunk))

    progress_bar.close()
