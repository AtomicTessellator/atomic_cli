# cli/core/client.py
import logging
import sys
from functools import wraps

from typing import Optional, Dict, Any, Iterator, List, Union
import httpx
from rich.console import Console

from .config import Config

console = Console(stderr=True)
logger = logging.getLogger(__name__)


def handle_connection_errors(func):
    @wraps(func)
    def wrapper(*args, **kwargs):
        try:
            return func(*args, **kwargs)
        except (httpx.ConnectError, httpx.ReadError, httpx.ReadTimeout) as e:
            console.print("[red]Connection error: Unable to communicate with the server.[/red]")
            console.print("[white]Make sure the server is running and accessible.[/white]")
            logger.error(f"Connection error: {str(e)}")
            sys.exit(1)
    return wrapper


class APIClient:
    def __init__(self, base_url: str = "http://localhost:5004"):
        self.base_url = base_url
        self.params = {'limit': 20}
        self.client = httpx.Client(base_url=base_url, timeout=30.0)
        self._token: Optional[str] = None

    def set_auth(self, username: str, password: str):
        """Set basic auth credentials"""
        self.client.auth = (username, password)

    def set_token(self, token: str):
        """Set bearer token"""
        self._token = token
        self.client.headers["Authorization"] = f"Token {token}"

    def _handle_response(self, response: httpx.Response) -> Union[List, Dict[str, Any]]:
        """Handle API response and common status codes"""
        try:
            response.raise_for_status()
            if response.request.method in ('DELETE', 'HEAD', 'TRACE'):
                # no JSON returned
                return
            return response.json()
        except httpx.ConnectError as e:
            console.print("[red]Unable to connect to the server. Please check if the server is running and accessible.[/red]")
            logger.debug(f"Connection error: {str(e)}")
            sys.exit(1)
        except httpx.ReadError as e:
            console.print("[red]Connection lost while communicating with the server.[/red]")
            logger.debug(f"Read error: {str(e)}")
            sys.exit(1)
        except httpx.HTTPStatusError as e:
            error_data = {}
            try:
                error_data = e.response.json()
            except ValueError:
                pass
            
            # TODO: standardize error message responses
            if isinstance(error_data, dict):
                error_message = error_data.get('errors', e.response.text)

            if e.response.status_code == 400:
                console.print("[red]Invalid request. Please check your input.[/red]")
                if isinstance(error_message, dict):
                    for field, errors in error_message.items():
                        console.print(f"[red]  {field}: {', '.join(errors)}[/red]")
            elif e.response.status_code == 401:
                console.print("[red]Authentication failed. Please check your credentials or log in again.[/red]")
            elif e.response.status_code == 403:
                console.print("[red]Permission denied. You don't have access to this resource.[/red]")
                console.print(f"[white]{error_data}")
            elif e.response.status_code == 404:
                console.print("[red]Resource not found. Please check the ID or path.[/red]")
            elif e.response.status_code == 429:
                console.print("[red]Too many requests. Please try again later.[/red]")
            elif e.response.status_code >= 500:
                console.print("[red]Server error. Please try again later or contact support.[/red]")

            logger.debug(f"Server error response: {error_data}")
            sys.exit(1)

    @handle_connection_errors
    def get(self, path: str, params: Optional[Dict] = None) -> Dict[str, Any]:
        """Make GET request"""
        params = {**self.params, **params} if params else self.params
        response = self.client.get(path, params=params)
        return self._handle_response(response)

    @handle_connection_errors
    def post(self, path: str, data: Dict) -> Dict[str, Any]:
        """Make POST request"""
        response = self.client.post(path, json=data)
        return self._handle_response(response)

    @handle_connection_errors
    def put(self, path: str, data: Dict) -> Dict[str, Any]:
        """Make PUT request"""
        response = self.client.put(path, json=data)
        return self._handle_response(response)
    
    @handle_connection_errors
    def delete(self, path: str) -> None:
        """Make DELETE request"""
        response = self.client.delete(path)
        self._handle_response(response)

    @handle_connection_errors
    def paginate(self, path: str, params: Optional[Dict] = None) -> Iterator[Dict[str, Any]]:
        """Handle paginated responses"""
        params = params or {}
        while path:
            response = self.get(path, params)
            if isinstance(response, list):
                raise RuntimeError("Not supported. Please contact support with this error.")
            elif isinstance(response, dict):
                for result in response.get('results', []):
                    yield result
                path = response.get('next')
                # If there's a next page, convert full URL back to path
                if path:
                    path = path.replace(self.base_url, '')
                params = {}  # Clear params for subsequent requests
            else:
                break

    def get_all(self, path: str, params: Optional[Dict] = None) -> List[Dict[str, Any]]:
        """Get all results from a paginated endpoint"""
        return list(self.paginate(path, params))


def get_client() -> APIClient:
    """Get a configured API client instance"""
    config = Config()
    client = APIClient()
    
    if config.token:
        client.set_token(config.token)
    else:
        config.ensure_auth()
        client.set_auth(config.username, config.password)
        
    return client
