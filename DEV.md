# Dev Notes

- Avoid enabling rankings polling until AF3 live-refresh is re-implemented. Previously we had an 8 s interval (`scheduleRankingsPolling` in `webapp/static/js/app.js`) that spammed the server and starved other async work when the websocket bridge was down.
- When debugging cluster interactions locally, temporarily disable `scheduleRankingsPolling()` (comment invocation or bump interval).

