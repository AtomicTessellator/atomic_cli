import logging

import pynvml

# Having this hard coded list really sucks, but I don't know of a better way to do it
# pycuda library adds more than 1GB to the size of the continer
# and deviceQuery is statically linked to the CUDA library and fails on our setup

# Mapping of GPU models to their compute capabilities
compute_capability_map = {
    "NVIDIA GeForce RTX 3060": "86",
    "Tesla P40": "61",
}


def gpu_compute_capacity():

    try:
        pynvml.nvmlInit()

        device_count = pynvml.nvmlDeviceGetCount()
        if device_count == 0:
            logging.info("No CUDA-capable devices found")
            return None

        for i in range(device_count):
            handle = pynvml.nvmlDeviceGetHandleByIndex(i)
            device_name = pynvml.nvmlDeviceGetName(handle).decode("utf-8")

            compute_capability = compute_capability_map.get(device_name, "Unknown")
            if compute_capability != "Unknown":
                return compute_capability
            else:
                logging.error(f"Compute capability for {device_name} is not in the map")
                return None

    except Exception as e:
        logging.error(f"Failed to extract compute capability: {e}")
        return None
    finally:
        pynvml.nvmlShutdown()
