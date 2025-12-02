import urllib.request, os, json

MAPPING_FILE = "downloaded_urls.json"

# === manage downloads with URL mapping ===
def download_input_genomes(urls):
    # Load mapping if exists
    if os.path.exists(MAPPING_FILE):
        with open(MAPPING_FILE, "r") as f:
            url_map = json.load(f)
    else:
        url_map = {}

    # Determine current max input number
    existing_files = [f for f in os.listdir() if f.startswith("input") and f.endswith(".fna.gz")]
    max_num = 0
    for f in existing_files:
        try:
            n = int(f.replace("input", "").replace(".fna.gz", ""))
            if n > max_num:
                max_num = n
        except:
            continue

    local_files = []
    next_num = max_num + 1

    for url in urls:
        if not url.endswith(".fna.gz"):
            raise ValueError(f"URL must end with .fna.gz.")

        if url in url_map:
            local_file = url_map[url]
            print("The URL is already saved and downloaded!")
            if not os.path.exists(local_file):
                # File missing even though URL is mapped
                print(f"{local_file} mapped to URL but missing locally. Downloading again...")
                urllib.request.urlretrieve(url, local_file)
                print("Download complete!")
        else:
            # Assign new input{i}.fna.gz
            local_file = f"input{next_num}.fna.gz"
            print(f"Downloading {url} -> {local_file}")
            urllib.request.urlretrieve(url, local_file)
            print("Download complete!")
            url_map[url] = local_file
            next_num += 1

        local_files.append(local_file)

    # Save updated mapping
    with open(MAPPING_FILE, "w") as f:
        json.dump(url_map, f, indent=2)

    return local_files
