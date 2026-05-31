"""Test that LokiHandler properly flushes all logs when close() is called"""

import logging
import time
from unittest.mock import patch, MagicMock

from atomict.infra.logging import LokiHandler


def test_loki_handler_flushes_batch_on_close():
    """
    Regression test: Ensures that logs in the worker's batch are flushed when close() is called.
    
    This test verifies the fix for a bug where logs sitting in the worker's internal batch
    were lost when close() was called, because:
    1. The worker thread wouldn't wake up from queue.get() until timeout
    2. The worker thread didn't flush its final batch before exiting
    """
    
    # Track all calls to requests.post
    post_calls = []
    
    def mock_post(*args, **kwargs):
        """Mock requests.post and record what was sent"""
        post_calls.append({
            'url': args[0] if args else kwargs.get('url'),
            'json': kwargs.get('json'),
            'timestamp': time.time()
        })
        response = MagicMock()
        response.raise_for_status = MagicMock()
        return response
    
    with patch('atomict.infra.logging.requests.post', side_effect=mock_post):
        # Create handler with:
        # - High batch_size (100) so it won't auto-flush
        # - Long flush_interval (60s) so it won't time-flush
        handler = LokiHandler(
            url="http://localhost:3100",
            task_id="test-task",
            batch_size=100,
            flush_interval=60.0
        )
        handler.setFormatter(logging.Formatter('%(message)s'))
        
        # Send 5 log messages (well below batch_size of 100)
        for i in range(5):
            record = logging.LogRecord(
                name="test",
                level=logging.INFO,
                pathname="",
                lineno=0,
                msg=f"Test message {i}",
                args=(),
                exc_info=None
            )
            handler.emit(record)
        
        # Give the worker thread time to pick up the logs into its batch
        time.sleep(0.1)
        
        # Close the handler - this should flush all logs
        handler.close()
        
        # Count how many log messages were actually sent to Loki
        total_logs_sent = 0
        for call in post_calls:
            if call['json'] and 'streams' in call['json']:
                for stream in call['json']['streams']:
                    total_logs_sent += len(stream['values'])
        
        # All 5 logs should have been sent
        assert total_logs_sent == 5, \
            f"Expected all 5 logs to be sent, but only {total_logs_sent} were sent. " \
            f"Logs may have been lost during close()!"


def test_loki_handler_flushes_queued_logs_on_close():
    """
    Test that logs still in the queue (not yet picked up by worker) are also flushed on close.
    """
    
    post_calls = []
    
    def mock_post(*args, **kwargs):
        """Mock requests.post and record what was sent"""
        post_calls.append({
            'url': args[0] if args else kwargs.get('url'),
            'json': kwargs.get('json'),
        })
        response = MagicMock()
        response.raise_for_status = MagicMock()
        return response
    
    with patch('atomict.infra.logging.requests.post', side_effect=mock_post):
        handler = LokiHandler(
            url="http://localhost:3100",
            task_id="test-task",
            batch_size=100,
            flush_interval=60.0
        )
        handler.setFormatter(logging.Formatter('%(message)s'))
        
        # Send 3 log messages
        for i in range(3):
            record = logging.LogRecord(
                name="test",
                level=logging.INFO,
                pathname="",
                lineno=0,
                msg=f"Test message {i}",
                args=(),
                exc_info=None
            )
            handler.emit(record)
        
        # Close immediately without waiting - some logs might still be in queue
        handler.close()
        
        # Count logs sent
        total_logs_sent = 0
        for call in post_calls:
            if call['json'] and 'streams' in call['json']:
                for stream in call['json']['streams']:
                    total_logs_sent += len(stream['values'])
        
        # All logs should be sent, whether they were in the queue or batch
        assert total_logs_sent == 3, \
            f"Expected all 3 logs to be sent, but only {total_logs_sent} were sent."


def test_loki_handler_respects_task_id_label():
    """Test that task_id is properly included in log labels"""
    
    post_calls = []
    
    def mock_post(*args, **kwargs):
        post_calls.append(kwargs.get('json'))
        response = MagicMock()
        response.raise_for_status = MagicMock()
        return response
    
    with patch('atomict.infra.logging.requests.post', side_effect=mock_post):
        handler = LokiHandler(
            url="http://localhost:3100",
            task_id="my-task-123",
            batch_size=1,  # Flush immediately
            flush_interval=0.1
        )
        handler.setFormatter(logging.Formatter('%(message)s'))
        
        record = logging.LogRecord(
            name="test",
            level=logging.INFO,
            pathname="",
            lineno=0,
            msg="Test message",
            args=(),
            exc_info=None
        )
        handler.emit(record)
        
        # Wait for flush
        time.sleep(0.2)
        handler.close()
        
        # Check that task_id label was included
        assert len(post_calls) > 0, "No POST calls made"
        payload = post_calls[0]
        assert 'streams' in payload
        assert len(payload['streams']) > 0
        
        stream = payload['streams'][0]
        assert 'stream' in stream
        assert stream['stream'].get('task_id') == 'my-task-123', \
            f"Expected task_id label to be 'my-task-123', got {stream['stream'].get('task_id')}"


if __name__ == "__main__":
    print("Running LokiHandler close() tests...\n")
    
    try:
        test_loki_handler_flushes_batch_on_close()
        print("✅ test_loki_handler_flushes_batch_on_close PASSED")
    except AssertionError as e:
        print(f"❌ test_loki_handler_flushes_batch_on_close FAILED: {e}")
    
    try:
        test_loki_handler_flushes_queued_logs_on_close()
        print("✅ test_loki_handler_flushes_queued_logs_on_close PASSED")
    except AssertionError as e:
        print(f"❌ test_loki_handler_flushes_queued_logs_on_close FAILED: {e}")
    
    try:
        test_loki_handler_respects_task_id_label()
        print("✅ test_loki_handler_respects_task_id_label PASSED")
    except AssertionError as e:
        print(f"❌ test_loki_handler_respects_task_id_label FAILED: {e}")
    
    print("\nAll tests completed!")
