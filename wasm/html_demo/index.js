'use strict';
import 'maplibre-gl/dist/maplibre-gl.css';
import * as maplibregl from 'maplibre-gl';
import { Protocol } from 'pmtiles';
import { Lrs, LrmScaleMeasure, set_panic_hook } from '../pkg/liblrs_wasm';
import * as turf from '@turf/helpers';
import Bbox from '@turf/bbox'
import Alpine from 'alpinejs'

// For the rust bindings: this allows us to have nice error messages
set_panic_hook()

async function file_selected(el) {
    const [file] = el.target.files;
    const data = await file.arrayBuffer()
    const lrs = await Lrs.load(new Uint8Array(data));

    const curves_features = []
    const anchors_features = []
    for (let i = 0; i < lrs.lrm_len(); i++) {
        const anchors = lrs.get_anchors(i);
        const lrm_id = lrs.get_lrm_scale_id(i);
        const geom = lrs.get_lrm_geom(i);
        const feature = turf.lineString(geom.map(p => [p.x, p.y]), { id: lrm_id, anchors: anchors.map((a) => a.name) }, {
            id: i,
        });
        feature.bbox = Bbox(feature)
        curves_features.push(feature);

        for (const anchor of anchors) {
            const properties = { id: anchor.name, lrm_id, curve: anchor.curve_position, scale: anchor.scale_position }
            anchors_features.push(turf.point([anchor.position.x, anchor.position.y], properties, { id: anchors_features.length }))
        }
    }

    map.addSource('lrms', {
        'type': 'geojson',
        'data': turf.featureCollection(curves_features),
    });

    map.addSource('anchors', {
        'type': 'geojson',
        'data': turf.featureCollection(anchors_features),
    })

    map.addSource('pr', {
        'type': 'geojson',
        'data': turf.featureCollection([])
    })

    map.addSource('range', {
        type: 'geojson',
        data: turf.featureCollection([])
    })

    map.addLayer({
        'id': 'lrms',
        'type': 'line',
        'source': 'lrms',
        'paint': {
            'line-color': '#888',
            'line-width': 2
        }
    });

    map.addLayer({
        'id': 'lrms-hover',
        'type': 'line',
        'source': 'lrms',
        'paint': {
            'line-color': 'red',
            'line-width': 3,
            'line-opacity': [
                'case',
                ['boolean', ['feature-state', 'hover'], false],
                1,
                0
            ]
        }
    });

    map.addLayer({
        id: 'anchors-unselected',
        type: 'circle',
        source: 'anchors',
        paint: {
            'circle-radius': 5,
            'circle-color': '#eee',
            'circle-opacity': ['case',
                ['boolean', ['feature-state', 'selected'], true],
                0,
                1,
            ]
        }
    })

    map.addLayer({
        id: 'anchors',
        type: 'circle',
        source: 'anchors',
        paint: {
            'circle-radius': 5,
            'circle-color': 'blue',
            'circle-opacity': ['case',
                ['boolean', ['feature-state', 'selected'], true],
                1,
                0,
            ]
        }
    })

    map.addLayer({
        id: 'anchors-labels',
        type: 'symbol',
        source: 'anchors',
        layout: {
            //'text-field': ['concat', ['get', 'id'], ['literal', '-'], ['get', 'lrm_id'], ['literal', ' curve_pos:'], ['get', 'curve'], ['literal', ' scale:'], ['get', 'scale']],
            'text-field': ['get', 'id'],
            'text-offset': [1, 0],
        },
        paint: {
            'text-halo-color': 'white',
            'text-halo-width': 2,
        }
    })

    map.addLayer({
        id: 'pr-outline',
        type: 'circle',
        source: 'pr',
        paint: {
            'circle-radius': 10,
            'circle-color': 'white',
        }
    })

    map.addLayer({
        id: 'pr',
        type: 'circle',
        source: 'pr',
        paint: {
            'circle-radius': 6,
            'circle-color': 'red',
        }
    })

    map.addLayer({
        'id': 'range',
        'type': 'line',
        'source': 'range',
        'paint': {
            'line-color': 'yellow',
            'line-width': 2
        }
    });

    return {
        features: curves_features,
        filename: file.name,
        filesize: data.byteLength,
        anchors_features,
        lrs,
    }
}

let protocol = new Protocol();
maplibregl.addProtocol('pmtiles', protocol.tile);

map = new maplibregl.Map({
    container: 'map', // container id
    style: process.env.MAPLIBRE_STYLE,
    center: [2.3469, 46.8589], // starting position [lng, lat]
    zoom: 5, // starting zoom,
});


window.Alpine = Alpine
Alpine.store('lrms', {
    status: 'waiting',
    error_text: null,
    lrms: [],
    lrs: null,
    filesize: null,
    filename: null,
    selectedFeature: null,
    pkStart: '',
    pkEnd: '',
    pkStartPoint: null,
    pkEndPoint: null,
    startMeasure: null,
    endMeasure: null,
    anchors: [],
    anchors_features: [],
    filter: '',

    async load(el) {
        try {
            this.status = 'loading'
            const lrs = await file_selected(el)
            this.lrms = lrs.features
            this.filename = lrs.filename
            this.filesize = lrs.filesize
            this.status = 'loaded'
            this.lrs = lrs.lrs
            this.anchors_features = lrs.anchors_features
        } catch (e) {
            this.status = 'error'
            this.error_text = e
        }

    },
    get waiting() { return this.status == 'waiting' },
    get loaded() { return this.status == 'loaded' },
    get error() { return this.status == 'error' },
    get loading() { return this.status == 'loading' },
    details(id) {
        this.unselect()
        this.selectedFeature = this.lrms[id];

        map.setFeatureState(
            { source: 'lrms', id },
            { hover: true }
        )
        map.fitBounds(this.lrms[id].bbox, { padding: 40 })
        const lrm_id = this.lrms[id].properties.id


        for (const anchor of this.anchors_features) {
            if (anchor.properties.lrm_id !== lrm_id) {
                map.setFeatureState({ source: 'anchors', id: anchor.id }, { selected: false })
            }
        }
    },
    unselect() {
        if (this.selectedFeature) {
            map.setFeatureState(
                { source: 'lrms', id: this.selectedFeature.id },
                { hover: false }
            )
            this.selectedFeature = null
            for (const anchor of this.anchors_features) {
                map.setFeatureState({ source: 'anchors', id: anchor.id }, { selected: true })
            }
        }
    },
    startPkChange({ target }) {
        const re = /([0-9]+)\+([0-9]+)/;
        if (re.test(target.value)) {
            const [all, anchor_name, scale_offset] = target.value.match(re)
            console.log('anchorname', anchor_name, 'scale offset', scale_offset)
            this.startMeasure = new LrmScaleMeasure(anchor_name, scale_offset);
            const point = this.lrs.resolve(this.selectedFeature.id, this.startMeasure)
            this.pkStartPoint = turf.point([point.x, point.y]);
            this.handlePks()
        } else {
            this.pkStartPoint = null;
            this.startMeasure = null;
        }
    },
    endPkChange({ target }) {
        const re = /([0-9]+)\+([0-9]+)/;
        if (re.test(target.value)) {
            const [all, anchor_name, scale_offset] = target.value.match(re)
            this.endMeasure = new LrmScaleMeasure(anchor_name, scale_offset);
            const point = this.lrs.resolve(this.selectedFeature.id, this.endMeasure)
            this.pkEndPoint = turf.point([point.x, point.y])
            this.handlePks()
        } else {
            this.pkEndPoint = null;
            this.endMeasure = null;
        }
    },
    handlePks() {
        const points = [this.pkStartPoint, this.pkEndPoint].filter(p => p !== null)
        const geojson = turf.featureCollection(points);
        map.getSource('pr').setData(geojson);
        if (points.length === 1) {
            map.flyTo({ center: points[0].geometry.coordinates, zoom: 15 })
        } else {
            map.fitBounds(Bbox(geojson), { padding: 30 })
            const range = this.lrs.resolve_range(this.selectedFeature.id, this.startMeasure, this.endMeasure)
            const feature = turf.lineString(range.map(p => [p.x, p.y]));
            map.getSource('range').setData(turf.featureCollection([feature]))
        }
    }
})
Alpine.start()
