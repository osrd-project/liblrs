const path = require('path');
const HtmlWebpackPlugin = require('html-webpack-plugin');
const webpack = require('webpack');
const WasmPackPlugin = require("@wasm-tool/wasm-pack-plugin");

module.exports = {
    entry: './html_demo/index.js',
    output: {
        path: path.resolve(__dirname, 'dist'),
        filename: 'index.js',
    },
    plugins: [
        new HtmlWebpackPlugin({ template: path.resolve(__dirname, '/html_demo/index.html') }),
        new WasmPackPlugin({
            crateDirectory: path.resolve(__dirname, ".."),
            forceMode: "production", // We want to force to build in release
            outDir: path.resolve(__dirname, 'pkg')
        }),
        // Have this example work in Edge which doesn't ship `TextEncoder` or
        // `TextDecoder` at this time.
        new webpack.ProvidePlugin({
            TextDecoder: ['text-encoding', 'TextDecoder'],
            TextEncoder: ['text-encoding', 'TextEncoder']
        }),
        new webpack.EnvironmentPlugin({
            MAPLIBRE_STYLE: 'https://demotiles.maplibre.org/style.json',
        })
    ],
    mode: 'development',
    experiments: {
        asyncWebAssembly: true
    },
    module: {
        rules: [
            {
                test: /\.(css)$/,
                use: ['style-loader', 'css-loader'],
            },
        ],
    },
};
