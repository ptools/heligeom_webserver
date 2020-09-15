#! /usr/bin/env python

from heligeom import create_app

if __name__ == "__main__":
    app = create_app('config.DevConfig')
    app.run(debug=True,host="0.0.0.0")
