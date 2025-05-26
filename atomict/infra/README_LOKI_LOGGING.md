# Loki Logging Integration

This module provides automatic logging to both stdout and Loki (when configured).

## Configuration

Set the `AT_LOGGING_ENDPOINT` environment variable to your Loki server URL:

```bash
export AT_LOGGING_ENDPOINT=http://spitfire:12345
```

## Usage

```python
from atomict.infra.logging import config_loggers
import logging

# Configure logging with optional task_id for filtering
config_loggers(prefix="[MyApp]", task_id="task-123")

# Use standard Python logging
logger = logging.getLogger(__name__)
logger.info("This will be logged to both stdout and Loki")
```

## Features

- **Dual Output**: Logs are always sent to stdout, and optionally to Loki if `AT_LOGGING_ENDPOINT` is set
- **Task ID Labeling**: Pass a `task_id` to `config_loggers()` to label all logs for easy filtering
- **Automatic Batching**: Logs are batched and sent to Loki every 5 seconds or when 100 logs accumulate
- **Non-blocking**: Loki logging happens in a background thread and won't block your application
- **Graceful Degradation**: If Loki is unavailable, logs still go to stdout and errors are printed to stderr

## Querying Logs in Loki

To query logs for a specific task:

```bash
curl -G http://spitfire:12345/loki/api/v1/query_range \
     -H "X-Scope-OrgID: admin" \
     --data-urlencode 'query={task_id="task-123"}' \
     --data-urlencode 'start='$(date --date="5 minutes ago" +%s)'000000000' \
     --data-urlencode 'end='$(date +%s)'000000000'
```

Other useful queries:
- All Python logs: `{job="python"}`
- Logs by level: `{job="python", level="error"}`
- Logs by logger name: `{job="python", logger="atomict.api"}`
- Combined filters: `{task_id="task-123", level=~"error|warning"}`

## Labels

Each log entry is labeled with:
- `job`: Always "python"
- `level`: Log level (debug, info, warning, error, critical)
- `logger`: Logger name (e.g., "atomict.api")
- `task_id`: Optional task identifier (if provided to `config_loggers`) 