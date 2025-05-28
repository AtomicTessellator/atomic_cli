import os
import time
import json
import logging
import threading
from queue import Queue
from logging.config import dictConfig
from typing import Optional
import requests


class LokiHandler(logging.Handler):
    """Custom logging handler that sends logs to Loki using requests"""
    
    def __init__(self, url: str, task_id: Optional[str] = None, batch_size: int = 100, flush_interval: float = 5.0, auth: Optional[tuple[str, str]] = None):
        super().__init__()
        self.url = url
        self.task_id = task_id
        self.batch_size = batch_size
        self.flush_interval = flush_interval
        self.queue = Queue()
        self.stop_event = threading.Event()
        
        if auth is None:
            # Try and get auth from environment variables
            username = os.environ.get("AT_LOGGING_USERNAME")
            password = os.environ.get("AT_LOGGING_PASSWORD")
            if username and password:
                auth = (username, password)

        # _auth is set once during init and never modified, making it safe to read from any thread
        self._auth = auth
        # Start background thread for batching and sending logs
        self.worker_thread = threading.Thread(target=self._worker, daemon=True)
        self.worker_thread.start()
    
    def emit(self, record):
        """Queue log record for sending to Loki"""
        try:
            msg = self.format(record)
            
            # Create Loki-compatible log entry
            timestamp = str(int(record.created * 1e9))  # nanoseconds
            
            # Build labels
            labels = {
                "job": "python",
                "level": record.levelname.lower(),
                "logger": record.name,
            }
            
            if self.task_id:
                labels["task_id"] = self.task_id
            
            # Add to queue
            entry = {
                "stream": labels,
                "values": [[timestamp, msg]]
            }
            self.queue.put(entry)
            
        except Exception as e:
            print(f"[LokiHandler.emit] ERROR: {e}")
            import traceback
            traceback.print_exc()
            self.handleError(record)
    
    def _worker(self):
        """Background worker that batches and sends logs to Loki"""
        
        try:
            batch = []
            last_flush = time.time()
            
            while not self.stop_event.is_set():
                try:
                    # Try to get items from queue with timeout
                    timeout = max(0.1, self.flush_interval - (time.time() - last_flush))
                    
                    try:
                        item = self.queue.get(timeout=timeout)
                        batch.append(item)
                    except:
                        pass  # Queue.get timed out
                    
                    # Check if we should flush
                    time_since_flush = time.time() - last_flush
                    should_flush = (
                        len(batch) >= self.batch_size or 
                        time_since_flush >= self.flush_interval
                    )
                    
                    if should_flush:
                        if batch:
                            self._send_batch(batch)
                            batch = []
                        else:
                            pass
                        last_flush = time.time()  # Update last_flush even if batch was empty
                        
                except Exception as e:
                    # Log to stderr to avoid recursion
                    import sys
                    print(f"[LokiHandler._worker] ERROR in worker: {e}", file=sys.stderr)
                    import traceback
                    traceback.print_exc()
            
        except Exception as e:
            import sys
            print(f"[LokiHandler._worker] FATAL ERROR in worker thread: {e}", file=sys.stderr)
            import traceback
            traceback.print_exc()
        
    def _send_batch(self, batch):
        """Send a batch of logs to Loki"""
        try:
            # Merge streams with same labels
            streams_dict = {}
            for item in batch:
                stream_key = json.dumps(item["stream"], sort_keys=True)
                if stream_key not in streams_dict:
                    streams_dict[stream_key] = {
                        "stream": item["stream"],
                        "values": []
                    }
                streams_dict[stream_key]["values"].extend(item["values"])
            
            # Prepare payload
            payload = {
                "streams": list(streams_dict.values())
            }
            
            # Send to Loki
            headers = {
                "Content-Type": "application/json",
                "X-Scope-OrgID": "admin"
            }
            
            response = requests.post(
                f"{self.url}/loki/api/v1/push",
                json=payload,
                headers=headers,
                timeout=5,
                auth=self._auth
            )
            response.raise_for_status()
            
        except Exception as e:
            # Log to stderr to avoid recursion
            import sys
            print(f"[LokiHandler._send_batch] ERROR sending to Loki: {e}", file=sys.stderr)
    
    def close(self):
        """Flush remaining logs and stop the worker thread"""
        
        self.stop_event.set()
        # Flush any remaining logs
        remaining = []
        while not self.queue.empty():
            try:
                remaining.append(self.queue.get_nowait())
            except:
                break
        if remaining:
            self._send_batch(remaining)
        
        # Wait for worker to finish
        self.worker_thread.join(timeout=1)
        super().close()


def config_loggers(
    prefix: str = '',
    task_id: Optional[str] = None,
    batch_size: int = 10,
    flush_interval: float = 2.0,
    *args,
    **kwargs
):
    """Configure loggers with console and optional Loki output
    
    Args:
        prefix: Prefix to add to log messages
        task_id: Task ID to include in Loki labels for filtering
    """
    # Base configuration with console handler
    logging_config = dict(
        version=1,
        formatters={
            "verbose": {"format": f"%(levelname)s %(asctime)s {prefix} %(message)s"}
        },
        handlers={
            "console": {"class": "logging.StreamHandler", "formatter": "verbose"}
        },
        root={"handlers": ["console"], "level": "INFO"},
        loggers={
            "atomict.api": {
                "handlers": ["console"],
                "level": "INFO",
                "propagate": False
            }
        }
    )
    
    # Apply the configuration FIRST
    dictConfig(logging_config)
    
    # Add Loki handler if endpoint is configured AFTER dictConfig
    loki_endpoint = os.environ.get("AT_LOGGING_ENDPOINT")
    
    if loki_endpoint:
        # Create and configure Loki handler
        loki_handler = LokiHandler(
            url=loki_endpoint,
            task_id=task_id,
            batch_size=batch_size,
            flush_interval=flush_interval,
        )
        loki_handler.setFormatter(logging.Formatter(f"%(levelname)s %(asctime)s {prefix} %(message)s"))
        
        # Add to root logger
        logging.getLogger().addHandler(loki_handler)
        
        # Also add to specific loggers
        for logger_name in ["atomict.api"]:
            logger = logging.getLogger(logger_name)
            logger.addHandler(loki_handler)
