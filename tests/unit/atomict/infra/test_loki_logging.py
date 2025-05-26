#!/usr/bin/env python3
"""
Test script for Loki logging integration

Usage:
    export AT_LOGGING_ENDPOINT=http://spitfire:12345
    python test_loki_logging.py
    
Then query Loki:
    curl -G http://spitfire:12345/loki/api/v1/query_range \
         -H "X-Scope-OrgID: admin" \
         --data-urlencode 'query={task_id="test-task-123"}' \
         --data-urlencode 'start='$(date --date="5 minutes ago" +%s)'000000000' \
         --data-urlencode 'end='$(date +%s)'000000000'
"""

import os
import logging
import time
from atomict.infra.logging import config_loggers

# Set up logging with a task_id
task_id = "test-task-123"
config_loggers(prefix="[TEST]", task_id=task_id)

# Get logger
logger = logging.getLogger(__name__)

# Test logging at different levels
logger.info(f"Starting test task with ID: {task_id}")
logger.debug("This is a debug message")
logger.warning("This is a warning message")

# Simulate some work
for i in range(5):
    logger.info(f"Processing item {i+1}/5")
    time.sleep(1)

logger.info("Task completed successfully")

# Test error logging
try:
    raise ValueError("Test error for demonstration")
except Exception as e:
    logger.error(f"Caught exception: {e}", exc_info=True)

# Give time for logs to be sent
time.sleep(2)
print("\nLogs have been sent to stdout and Loki (if AT_LOGGING_ENDPOINT is set)")
print(f"Query Loki with: {{task_id=\"{task_id}\"}}") 