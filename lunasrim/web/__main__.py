"""Start the LunaSRIM web application.

    python -m lunasrim.web [--host 127.0.0.1] [--port 8000]
"""
import argparse


def main(argv=None):
    p = argparse.ArgumentParser(prog='lunasrim.web', description='LunaSRIM web application')
    p.add_argument('--host', default='127.0.0.1')
    p.add_argument('--port', type=int, default=8000)
    a = p.parse_args(argv)
    try:
        import uvicorn
    except ImportError:
        raise SystemExit('the web application needs fastapi and uvicorn: pip install "lunasrim[web]"')
    print('LunaSRIM web application on http://%s:%d/' % (a.host, a.port))
    uvicorn.run('lunasrim.web.server:app', host=a.host, port=a.port)


if __name__ == '__main__':
    main()
